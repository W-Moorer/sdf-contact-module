// ============================================================================
// SdfOracle.cpp — Query engine implementation
// ============================================================================

#include "SdfOracle/SdfOracle.h"
#include "SdfOracle/SdfIO.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sdf {

namespace {
struct SparseSampleCache
{
    explicit SparseSampleCache(const SparseSdfData& source)
        : data(source)
    {
    }

    const SparseSdfBlock* blockForVoxel(int i, int j, int k)
    {
        if (!data.spec.inBounds(i, j, k) || data.blockDim <= 0)
        {
            return nullptr;
        }
        const int bi = i / data.blockDim;
        const int bj = j / data.blockDim;
        const int bk = k / data.blockDim;
        const uint64_t key = SparseSdfData::key(bi, bj, bk);
        for (int n = 0; n < count; ++n)
        {
            if (keys[static_cast<size_t>(n)] == key)
            {
                return blocks[static_cast<size_t>(n)];
            }
        }

        const auto it = data.blockLookup.find(key);
        const SparseSdfBlock* block =
            it == data.blockLookup.end() ? nullptr : &data.blocks[it->second];
        if (count < static_cast<int>(keys.size()))
        {
            keys[static_cast<size_t>(count)] = key;
            blocks[static_cast<size_t>(count)] = block;
            ++count;
        }
        return block;
    }

    bool sample(
        int i,
        int j,
        int k,
        double& phi,
        Eigen::Vector3d& normal,
        Eigen::Matrix3d* hessian = nullptr)
    {
        const SparseSdfBlock* block = blockForVoxel(i, j, k);
        return block != nullptr &&
            data.sampleFromBlock(*block, i, j, k, phi, normal, hessian);
    }

    const SparseSdfData& data;
    std::array<uint64_t, 8> keys{};
    std::array<const SparseSdfBlock*, 8> blocks{};
    int count{0};
};

SdfQueryResult withPrimaryCandidate(SdfQueryResult result)
{
    if (result.valid && result.candidates.empty())
    {
        result.candidates.push_back(SdfQueryCandidate{
            result.gap,
            result.gapNormal,
            result.witnessPoint,
            0});
    }
    return result;
}

SdfFastQueryResult fastFromQueryResult(const SdfQueryResult& result)
{
    SdfFastQueryResult fast;
    fast.gap = result.gap;
    fast.gapNormal = result.gapNormal;
    fast.witnessPoint = result.witnessPoint;
    fast.valid = result.valid;
    if (!result.valid)
    {
        return fast;
    }

    if (result.candidates.empty())
    {
        fast.candidateCount = 1;
        fast.candidates[0] = SdfQueryCandidate{
            result.gap,
            result.gapNormal,
            result.witnessPoint,
            0};
        return fast;
    }

    fast.candidateCount = std::min(
        static_cast<int>(result.candidates.size()),
        kMaxSdfQueryCandidates);
    for (int i = 0; i < fast.candidateCount; ++i)
    {
        fast.candidates[static_cast<size_t>(i)] =
            result.candidates[static_cast<size_t>(i)];
    }
    return fast;
}

SdfQueryResult queryResultFromFast(const SdfFastQueryResult& fast)
{
    SdfQueryResult result;
    result.gap = fast.gap;
    result.gapNormal = fast.gapNormal;
    result.witnessPoint = fast.witnessPoint;
    result.valid = fast.valid;
    if (!fast.valid)
    {
        return result;
    }

    const int count = std::max(0, std::min(fast.candidateCount, kMaxSdfQueryCandidates));
    result.candidates.reserve(static_cast<size_t>(std::max(1, count)));
    for (int i = 0; i < count; ++i)
    {
        result.candidates.push_back(fast.candidates[static_cast<size_t>(i)]);
    }
    if (result.candidates.empty())
    {
        result.candidates.push_back(SdfQueryCandidate{
            result.gap,
            result.gapNormal,
            result.witnessPoint,
            0});
    }
    return result;
}
} // namespace

// ---------- Construction ----------

SdfOracle::SdfOracle(SdfData data)
    : m_data(std::move(data))
{
}

SdfOracle::SdfOracle(SparseSdfData data)
    : m_sparseData(std::move(data))
    , m_sparse(true)
{
}

SdfOracle SdfOracle::loadFromFile(const std::string& path)
{
    if (isSparseSdf(path))
    {
        return SdfOracle(loadSparseSdf(path));
    }
    return SdfOracle(loadSdf(path));
}

// ---------- Public query API ----------

SdfQueryResult SdfOracle::query(
    const Eigen::Vector3d& point, QueryMode mode) const
{
    if (!isValid())
    {
        SdfQueryResult r;
        r.valid = false;
        return r;
    }

    if (m_sparse)
    {
        switch (mode)
        {
        case QueryMode::ContactAware:
            return m_sparseData.hasContactCandidates()
                ? queryResultFromFast(queryFast(point, mode))
                : withPrimaryCandidate(sparseTricubicQuery(point));
        case QueryMode::SecondOrder:
            return withPrimaryCandidate(sparseSecondOrderQuery(point));
        case QueryMode::FirstOrder:
            return withPrimaryCandidate(sparseFirstOrderQuery(point));
        case QueryMode::Trilinear:
            return withPrimaryCandidate(sparseTrilinearQuery(point));
        case QueryMode::Tricubic:
            return withPrimaryCandidate(sparseTricubicQuery(point));
        default:
            return withPrimaryCandidate(sparseSecondOrderQuery(point));
        }
    }

    switch (mode)
    {
    case QueryMode::ContactAware:
        return m_data.hasContactCandidates()
            ? queryResultFromFast(queryFast(point, mode))
            : withPrimaryCandidate(tricubicQuery(point));
    case QueryMode::SecondOrder:
        return withPrimaryCandidate(secondOrderQuery(point));
    case QueryMode::FirstOrder:
        return withPrimaryCandidate(firstOrderQuery(point));
    case QueryMode::Trilinear:
        return withPrimaryCandidate(trilinearQuery(point));
    case QueryMode::Tricubic:
        return withPrimaryCandidate(tricubicQuery(point));
    default:
        return withPrimaryCandidate(secondOrderQuery(point));
    }
}

SdfFastQueryResult SdfOracle::queryFast(
    const Eigen::Vector3d& point,
    QueryMode mode) const
{
    if (!isValid())
    {
        return {};
    }

    if (mode != QueryMode::ContactAware)
    {
        if (m_sparse)
        {
            switch (mode)
            {
            case QueryMode::SecondOrder:
                return fastFromQueryResult(sparseSecondOrderQuery(point));
            case QueryMode::FirstOrder:
                return fastFromQueryResult(sparseFirstOrderQuery(point));
            case QueryMode::Trilinear:
                return fastFromQueryResult(sparseTrilinearQuery(point));
            case QueryMode::Tricubic:
                return fastFromQueryResult(sparseTricubicQuery(point));
            default:
                return fastFromQueryResult(sparseSecondOrderQuery(point));
            }
        }

        switch (mode)
        {
        case QueryMode::SecondOrder:
            return fastFromQueryResult(secondOrderQuery(point));
        case QueryMode::FirstOrder:
            return fastFromQueryResult(firstOrderQuery(point));
        case QueryMode::Trilinear:
            return fastFromQueryResult(trilinearQuery(point));
        case QueryMode::Tricubic:
            return fastFromQueryResult(tricubicQuery(point));
        default:
            return fastFromQueryResult(secondOrderQuery(point));
        }
    }

    if (m_sparse)
    {
        if (!m_sparseData.hasContactCandidates())
        {
            return fastFromQueryResult(sparseTricubicQuery(point));
        }

        const auto& spec = m_sparseData.spec;
        int i, j, k;
        spec.indexFromWorld(point, i, j, k);
        SdfFastQueryResult result;
        const SparseSdfBlock* block = m_sparseData.blockForVoxel(i, j, k);
        if (block == nullptr)
        {
            result.valid = false;
            return result;
        }

        const int lx = i - block->index[0] * m_sparseData.blockDim;
        const int ly = j - block->index[1] * m_sparseData.blockDim;
        const int lz = k - block->index[2] * m_sparseData.blockDim;
        if (lx < 0 || lx >= m_sparseData.blockDim ||
            ly < 0 || ly >= m_sparseData.blockDim ||
            lz < 0 || lz >= m_sparseData.blockDim)
        {
            result.valid = false;
            return result;
        }

        const int local = block->localIndex(lx, ly, lz, m_sparseData.blockDim);
        const int count = std::min(
            static_cast<int>(block->contactCandidateCount[static_cast<size_t>(local)]),
            m_sparseData.contactCandidateStride);
        if (count <= 0)
        {
            return fastFromQueryResult(sparseFirstOrderQuery(point));
        }

        const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
        const Eigen::Vector3d dx = point - center;
        int selected = 0;
        double selectedAbsGap = std::numeric_limits<double>::max();
        result.candidateCount = std::min(count, kMaxSdfQueryCandidates);
        for (int candidate = 0; candidate < result.candidateCount; ++candidate)
        {
            const int idx = block->contactCandidateIndex(
                local,
                candidate,
                m_sparseData.contactCandidateStride);
            Eigen::Vector3d normal{
                static_cast<double>(block->contactCandidateNormal0[
                    static_cast<size_t>(3 * idx + 0)]),
                static_cast<double>(block->contactCandidateNormal0[
                    static_cast<size_t>(3 * idx + 1)]),
                static_cast<double>(block->contactCandidateNormal0[
                    static_cast<size_t>(3 * idx + 2)])};
            const double normalNorm = normal.norm();
            if (normalNorm > 1.0e-15)
            {
                normal /= normalNorm;
            }
            else
            {
                normal = Eigen::Vector3d::UnitX();
            }
            const double gap =
                static_cast<double>(block->contactCandidatePhi0[
                    static_cast<size_t>(idx)]) +
                normal.dot(dx);
            const Eigen::Vector3d witness = point - gap * normal;
            result.candidates[static_cast<size_t>(candidate)] =
                SdfQueryCandidate{
                    gap,
                    normal,
                    witness,
                    block->contactCandidateBranchId[static_cast<size_t>(idx)]};

            const double absGap = std::abs(gap);
            if (absGap < selectedAbsGap)
            {
                selectedAbsGap = absGap;
                selected = candidate;
            }
        }

        const SdfQueryCandidate& primary =
            result.candidates[static_cast<size_t>(selected)];
        result.gap = primary.gap;
        result.gapNormal = primary.gapNormal;
        result.witnessPoint = primary.witnessPoint;
        result.valid = true;
        return result;
    }

    if (!m_data.hasContactCandidates())
    {
        return fastFromQueryResult(tricubicQuery(point));
    }

    const auto& spec = m_data.spec;
    int i, j, k;
    spec.indexFromWorld(point, i, j, k);
    SdfFastQueryResult result;
    if (!spec.inBounds(i, j, k))
    {
        result.valid = false;
        return result;
    }

    const int voxel = spec.linearIndex(i, j, k);
    const int count = std::min(
        static_cast<int>(m_data.contactCandidateCount[static_cast<size_t>(voxel)]),
        m_data.contactCandidateStride);
    if (count <= 0)
    {
        return fastFromQueryResult(firstOrderQuery(point));
    }

    const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
    const Eigen::Vector3d dx = point - center;
    int selected = 0;
    double selectedAbsGap = std::numeric_limits<double>::max();
    result.candidateCount = std::min(count, kMaxSdfQueryCandidates);
    for (int candidate = 0; candidate < result.candidateCount; ++candidate)
    {
        const int idx = m_data.contactCandidateIndex(voxel, candidate);
        Eigen::Vector3d normal =
            m_data.contactCandidateNormal0[static_cast<size_t>(idx)];
        const double normalNorm = normal.norm();
        if (normalNorm > 1.0e-15)
        {
            normal /= normalNorm;
        }
        else
        {
            normal = Eigen::Vector3d::UnitX();
        }
        const double gap =
            m_data.contactCandidatePhi0[static_cast<size_t>(idx)] +
            normal.dot(dx);
        const Eigen::Vector3d witness = point - gap * normal;
        result.candidates[static_cast<size_t>(candidate)] =
            SdfQueryCandidate{
                gap,
                normal,
                witness,
                m_data.contactCandidateBranchId[static_cast<size_t>(idx)]};

        const double absGap = std::abs(gap);
        if (absGap < selectedAbsGap)
        {
            selectedAbsGap = absGap;
            selected = candidate;
        }
    }

    const SdfQueryCandidate& primary =
        result.candidates[static_cast<size_t>(selected)];
    result.gap = primary.gap;
    result.gapNormal = primary.gapNormal;
    result.witnessPoint = primary.witnessPoint;
    result.valid = true;
    return result;
}

std::vector<SdfQueryResult> SdfOracle::queryBatch(
    const std::vector<Eigen::Vector3d>& points, QueryMode mode) const
{
    std::vector<SdfQueryResult> results(points.size());
#pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(points.size()); ++i)
    {
        results[static_cast<size_t>(i)] = query(points[static_cast<size_t>(i)], mode);
    }
    return results;
}

// ---------- 2nd-order Taylor expansion ----------

SdfQueryResult SdfOracle::secondOrderQuery(
    const Eigen::Vector3d& p) const
{
    const auto& spec = m_data.spec;
    int i, j, k;
    spec.indexFromWorld(p, i, j, k);

    SdfQueryResult result;
    if (!spec.inBounds(i, j, k))
    {
        result.valid = false;
        return result;
    }

    // Voxel center
    Eigen::Vector3d vj = spec.centerFromIndex(i, j, k);
    Eigen::Vector3d dx = p - vj;

    // Precomputed data at this voxel
    double phi0 = m_data.phi(i, j, k);
    const Eigen::Vector3d& n0 = m_data.normal(i, j, k);
    const Eigen::Matrix3d& H0 = m_data.hessian(i, j, k);

    // 2nd-order Taylor expansion:
    //   φ(p) = φ₀ + n₀ᵀ·δx + ½·δxᵀ·H₀·δx
    double firstOrder = n0.dot(dx);
    double secondOrder = 0.5 * dx.transpose() * H0 * dx;
    result.gap = phi0 + firstOrder + secondOrder;

    // Updated normal: n(p) = normalize(n₀ + H₀·δx)
    Eigen::Vector3d updatedNormal = n0 + H0 * dx;
    double nLen = updatedNormal.norm();
    if (nLen > 1e-15)
        result.gapNormal = updatedNormal / nLen;
    else
        result.gapNormal = n0;

    // Witness point estimate: p - φ(p) * n(p)
    result.witnessPoint = p - result.gap * result.gapNormal;

    result.valid = true;
    return result;
}

// ---------- 1st-order Taylor expansion ----------

SdfQueryResult SdfOracle::firstOrderQuery(
    const Eigen::Vector3d& p) const
{
    const auto& spec = m_data.spec;
    int i, j, k;
    spec.indexFromWorld(p, i, j, k);

    SdfQueryResult result;
    if (!spec.inBounds(i, j, k))
    {
        result.valid = false;
        return result;
    }

    Eigen::Vector3d vj = spec.centerFromIndex(i, j, k);
    Eigen::Vector3d dx = p - vj;

    double phi0 = m_data.phi(i, j, k);
    const Eigen::Vector3d& n0 = m_data.normal(i, j, k);

    // 1st-order Taylor: φ(p) = φ₀ + n₀ᵀ·δx
    result.gap = phi0 + n0.dot(dx);
    result.gapNormal = n0;  // Normal fixed (no update)
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

// ---------- Trilinear interpolation ----------

SdfQueryResult SdfOracle::trilinearQuery(
    const Eigen::Vector3d& p) const
{
    const auto& spec = m_data.spec;

    // Compute fractional grid coordinates
    Eigen::Vector3d f = (p - spec.bmin).cwiseQuotient(spec.voxel)
        - Eigen::Vector3d(0.5, 0.5, 0.5);

    int i0 = static_cast<int>(std::floor(f.x()));
    int j0 = static_cast<int>(std::floor(f.y()));
    int k0 = static_cast<int>(std::floor(f.z()));

    double u = f.x() - i0;
    double v = f.y() - j0;
    double w = f.z() - k0;

    // Check all 8 corners are in bounds
    SdfQueryResult result;
    if (i0 < 0 || i0 + 1 >= spec.nx() ||
        j0 < 0 || j0 + 1 >= spec.ny() ||
        k0 < 0 || k0 + 1 >= spec.nz())
    {
        result.valid = false;
        return result;
    }

    // Trilinear interpolation of phi values
    double c000 = m_data.phi(i0,   j0,   k0);
    double c100 = m_data.phi(i0+1, j0,   k0);
    double c010 = m_data.phi(i0,   j0+1, k0);
    double c110 = m_data.phi(i0+1, j0+1, k0);
    double c001 = m_data.phi(i0,   j0,   k0+1);
    double c101 = m_data.phi(i0+1, j0,   k0+1);
    double c011 = m_data.phi(i0,   j0+1, k0+1);
    double c111 = m_data.phi(i0+1, j0+1, k0+1);

    double c00 = c000 * (1 - u) + c100 * u;
    double c10 = c010 * (1 - u) + c110 * u;
    double c01 = c001 * (1 - u) + c101 * u;
    double c11 = c011 * (1 - u) + c111 * u;

    double c0 = c00 * (1 - v) + c10 * v;
    double c1 = c01 * (1 - v) + c11 * v;

    result.gap = c0 * (1 - w) + c1 * w;

    // Gradient from finite differences of interpolated values
    double dfdx = ((c100 - c000) * (1-v) * (1-w) +
                   (c110 - c010) * v * (1-w) +
                   (c101 - c001) * (1-v) * w +
                   (c111 - c011) * v * w) / spec.voxel.x();
    double dfdy = ((c010 - c000) * (1-u) * (1-w) +
                   (c110 - c100) * u * (1-w) +
                   (c011 - c001) * (1-u) * w +
                   (c111 - c101) * u * w) / spec.voxel.y();
    double dfdz = ((c001 - c000) * (1-u) * (1-v) +
                   (c101 - c100) * u * (1-v) +
                   (c011 - c010) * (1-u) * v +
                   (c111 - c110) * u * v) / spec.voxel.z();

    Eigen::Vector3d grad(dfdx, dfdy, dfdz);
    double gLen = grad.norm();
    if (gLen > 1e-15)
        result.gapNormal = grad / gLen;
    else
        result.gapNormal = Eigen::Vector3d::UnitX();

    // Witness estimate
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

// ---------- Sparse block queries ----------

SdfQueryResult SdfOracle::sparseFirstOrderQuery(const Eigen::Vector3d& p) const
{
    const auto& spec = m_sparseData.spec;
    int i, j, k;
    spec.indexFromWorld(p, i, j, k);

    SdfQueryResult result;
    double phi0 = 0.0;
    Eigen::Vector3d n0;
    if (!m_sparseData.sample(i, j, k, phi0, n0, nullptr))
    {
        result.valid = false;
        return result;
    }

    const Eigen::Vector3d vj = spec.centerFromIndex(i, j, k);
    const Eigen::Vector3d dx = p - vj;
    result.gap = phi0 + n0.dot(dx);
    result.gapNormal = n0;
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

SdfQueryResult SdfOracle::sparseSecondOrderQuery(const Eigen::Vector3d& p) const
{
    if (!m_sparseData.hasHessian)
    {
        return sparseFirstOrderQuery(p);
    }

    const auto& spec = m_sparseData.spec;
    int i, j, k;
    spec.indexFromWorld(p, i, j, k);

    SdfQueryResult result;
    double phi0 = 0.0;
    Eigen::Vector3d n0;
    Eigen::Matrix3d H0;
    if (!m_sparseData.sample(i, j, k, phi0, n0, &H0))
    {
        result.valid = false;
        return result;
    }

    const Eigen::Vector3d vj = spec.centerFromIndex(i, j, k);
    const Eigen::Vector3d dx = p - vj;
    result.gap = phi0 + n0.dot(dx) + 0.5 * dx.transpose() * H0 * dx;
    const Eigen::Vector3d updatedNormal = n0 + H0 * dx;
    const double nLen = updatedNormal.norm();
    result.gapNormal = nLen > 1e-15 ? updatedNormal / nLen : n0;
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

SdfQueryResult SdfOracle::sparseTrilinearQuery(const Eigen::Vector3d& p) const
{
    const auto& spec = m_sparseData.spec;
    const Eigen::Vector3d f = (p - spec.bmin).cwiseQuotient(spec.voxel)
        - Eigen::Vector3d(0.5, 0.5, 0.5);

    const int i0 = static_cast<int>(std::floor(f.x()));
    const int j0 = static_cast<int>(std::floor(f.y()));
    const int k0 = static_cast<int>(std::floor(f.z()));
    const double u = f.x() - i0;
    const double v = f.y() - j0;
    const double w = f.z() - k0;

    SdfQueryResult result;
    if (i0 < 0 || i0 + 1 >= spec.nx() ||
        j0 < 0 || j0 + 1 >= spec.ny() ||
        k0 < 0 || k0 + 1 >= spec.nz())
    {
        result.valid = false;
        return result;
    }

    SparseSampleCache cache(m_sparseData);
    auto phiAt = [&](int i, int j, int k, double& phi) {
        Eigen::Vector3d n;
        return cache.sample(i, j, k, phi, n, nullptr);
    };

    double c000, c100, c010, c110, c001, c101, c011, c111;
    if (!phiAt(i0, j0, k0, c000) ||
        !phiAt(i0 + 1, j0, k0, c100) ||
        !phiAt(i0, j0 + 1, k0, c010) ||
        !phiAt(i0 + 1, j0 + 1, k0, c110) ||
        !phiAt(i0, j0, k0 + 1, c001) ||
        !phiAt(i0 + 1, j0, k0 + 1, c101) ||
        !phiAt(i0, j0 + 1, k0 + 1, c011) ||
        !phiAt(i0 + 1, j0 + 1, k0 + 1, c111))
    {
        result.valid = false;
        return result;
    }

    const double c00 = c000 * (1 - u) + c100 * u;
    const double c10 = c010 * (1 - u) + c110 * u;
    const double c01 = c001 * (1 - u) + c101 * u;
    const double c11 = c011 * (1 - u) + c111 * u;
    const double c0 = c00 * (1 - v) + c10 * v;
    const double c1 = c01 * (1 - v) + c11 * v;
    result.gap = c0 * (1 - w) + c1 * w;

    const double dfdx = ((c100 - c000) * (1-v) * (1-w) +
                         (c110 - c010) * v * (1-w) +
                         (c101 - c001) * (1-v) * w +
                         (c111 - c011) * v * w) / spec.voxel.x();
    const double dfdy = ((c010 - c000) * (1-u) * (1-w) +
                         (c110 - c100) * u * (1-w) +
                         (c011 - c001) * (1-u) * w +
                         (c111 - c101) * u * w) / spec.voxel.y();
    const double dfdz = ((c001 - c000) * (1-u) * (1-v) +
                         (c101 - c100) * u * (1-v) +
                         (c011 - c010) * (1-u) * v +
                         (c111 - c110) * u * v) / spec.voxel.z();

    const Eigen::Vector3d grad(dfdx, dfdy, dfdz);
    const double gLen = grad.norm();
    if (gLen > 1e-15)
        result.gapNormal = grad / gLen;
    else
        result.gapNormal = Eigen::Vector3d::UnitX();
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

SdfQueryResult SdfOracle::sparseTricubicQuery(const Eigen::Vector3d& p) const
{
    const auto& spec = m_sparseData.spec;
    const Eigen::Vector3d f = (p - spec.bmin).cwiseQuotient(spec.voxel)
        - Eigen::Vector3d(0.5, 0.5, 0.5);

    const int i0 = static_cast<int>(std::floor(f.x()));
    const int j0 = static_cast<int>(std::floor(f.y()));
    const int k0 = static_cast<int>(std::floor(f.z()));
    const double u = f.x() - i0;
    const double v = f.y() - j0;
    const double w = f.z() - k0;

    SdfQueryResult result;
    if (i0 < 0 || i0 + 1 >= spec.nx() ||
        j0 < 0 || j0 + 1 >= spec.ny() ||
        k0 < 0 || k0 + 1 >= spec.nz())
    {
        result.valid = false;
        return result;
    }

    double phi[2][2][2];
    double gx[2][2][2], gy[2][2][2], gz[2][2][2];
    const double vx = spec.voxel.x(), vy = spec.voxel.y(), vz = spec.voxel.z();
    SparseSampleCache cache(m_sparseData);

    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        Eigen::Vector3d n;
        if (!cache.sample(i0 + di, j0 + dj, k0 + dk, phi[di][dj][dk], n, nullptr))
        {
            result.valid = false;
            return result;
        }
        gx[di][dj][dk] = n.x() * vx;
        gy[di][dj][dk] = n.y() * vy;
        gz[di][dj][dk] = n.z() * vz;
    }

    auto h00 = [](double t) { return t * t * (2.0 * t - 3.0) + 1.0; };
    auto h10 = [](double t) { return t * (t * (t - 2.0) + 1.0); };
    auto h01 = [](double t) { return t * t * (3.0 - 2.0 * t); };
    auto h11 = [](double t) { return t * t * (t - 1.0); };

    const double a0 = h00(u), a1 = h10(u), a2 = h01(u), a3 = h11(u);
    const double b0 = h00(v), b1 = h10(v), b2 = h01(v), b3 = h11(v);
    const double c0 = h00(w), c1 = h10(w), c2 = h01(w), c3 = h11(w);

    double gap = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        const double p0 = (di == 0) ? a0 : a2;
        const double p1 = (di == 0) ? a1 : a3;
        const double q0 = (dj == 0) ? b0 : b2;
        const double q1 = (dj == 0) ? b1 : b3;
        const double r0 = (dk == 0) ? c0 : c2;
        const double r1 = (dk == 0) ? c1 : c3;
        gap += (p0 * q0 * r0) * phi[di][dj][dk];
        gap += (p1 * q0 * r0) * gx[di][dj][dk];
        gap += (p0 * q1 * r0) * gy[di][dj][dk];
        gap += (p0 * q0 * r1) * gz[di][dj][dk];
    }
    result.gap = gap;

    auto dh00 = [](double t) { return 6.0 * t * (t - 1.0); };
    auto dh10 = [](double t) { return t * (3.0 * t - 4.0) + 1.0; };
    auto dh01 = [](double t) { return 6.0 * t * (1.0 - t); };
    auto dh11 = [](double t) { return t * (3.0 * t - 2.0); };

    const double da0 = dh00(u), da1 = dh10(u), da2 = dh01(u), da3 = dh11(u);
    const double db0 = dh00(v), db1 = dh10(v), db2 = dh01(v), db3 = dh11(v);
    const double dc0 = dh00(w), dc1 = dh10(w), dc2 = dh01(w), dc3 = dh11(w);

    double dfdu = 0.0, dfdv = 0.0, dfdw = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        const double p0 = (di == 0) ? a0 : a2;
        const double p1 = (di == 0) ? a1 : a3;
        const double dp0 = (di == 0) ? da0 : da2;
        const double dp1 = (di == 0) ? da1 : da3;
        const double q0 = (dj == 0) ? b0 : b2;
        const double q1 = (dj == 0) ? b1 : b3;
        const double dq0 = (dj == 0) ? db0 : db2;
        const double dq1 = (dj == 0) ? db1 : db3;
        const double r0 = (dk == 0) ? c0 : c2;
        const double r1 = (dk == 0) ? c1 : c3;
        const double dr0 = (dk == 0) ? dc0 : dc2;
        const double dr1 = (dk == 0) ? dc1 : dc3;

        dfdu += (dp0 * q0 * r0) * phi[di][dj][dk];
        dfdu += (dp1 * q0 * r0) * gx[di][dj][dk];
        dfdu += (dp0 * q1 * r0) * gy[di][dj][dk];
        dfdu += (dp0 * q0 * r1) * gz[di][dj][dk];

        dfdv += (p0 * dq0 * r0) * phi[di][dj][dk];
        dfdv += (p1 * dq0 * r0) * gx[di][dj][dk];
        dfdv += (p0 * dq1 * r0) * gy[di][dj][dk];
        dfdv += (p0 * dq0 * r1) * gz[di][dj][dk];

        dfdw += (p0 * q0 * dr0) * phi[di][dj][dk];
        dfdw += (p1 * q0 * dr0) * gx[di][dj][dk];
        dfdw += (p0 * q1 * dr0) * gy[di][dj][dk];
        dfdw += (p0 * q0 * dr1) * gz[di][dj][dk];
    }

    const Eigen::Vector3d grad(dfdu / vx, dfdv / vy, dfdw / vz);
    const double gLen = grad.norm();
    if (gLen > 1e-15)
        result.gapNormal = grad / gLen;
    else
        result.gapNormal = Eigen::Vector3d::UnitX();
    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

SdfQueryResult SdfOracle::contactAwareQuery(const Eigen::Vector3d& p) const
{
    const auto& spec = m_data.spec;
    int i, j, k;
    spec.indexFromWorld(p, i, j, k);

    SdfQueryResult result;
    if (!spec.inBounds(i, j, k) || !m_data.hasContactCandidates())
    {
        result.valid = false;
        return result;
    }

    const int voxel = spec.linearIndex(i, j, k);
    const int count = static_cast<int>(m_data.contactCandidateCount[voxel]);
    if (count <= 0)
    {
        return firstOrderQuery(p);
    }

    const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
    const Eigen::Vector3d dx = p - center;

    int selected = 0;
    double selectedAbsGap = std::numeric_limits<double>::max();
    result.candidates.reserve(static_cast<size_t>(count));
    for (int candidate = 0; candidate < count; ++candidate)
    {
        const int idx = m_data.contactCandidateIndex(voxel, candidate);
        Eigen::Vector3d normal = m_data.contactCandidateNormal0[idx];
        const double normalNorm = normal.norm();
        if (normalNorm > 1.0e-15)
        {
            normal /= normalNorm;
        }
        else
        {
            normal = Eigen::Vector3d::UnitX();
        }
        const double gap = m_data.contactCandidatePhi0[idx] + normal.dot(dx);
        const Eigen::Vector3d witness = p - gap * normal;
        result.candidates.push_back(SdfQueryCandidate{
            gap,
            normal,
            witness,
            m_data.contactCandidateBranchId[static_cast<size_t>(idx)]});

        const double absGap = std::abs(gap);
        if (absGap < selectedAbsGap)
        {
            selectedAbsGap = absGap;
            selected = candidate;
        }
    }

    const SdfQueryCandidate& primary = result.candidates[static_cast<size_t>(selected)];
    result.gap = primary.gap;
    result.gapNormal = primary.gapNormal;
    result.witnessPoint = primary.witnessPoint;
    result.valid = true;
    return result;
}

SdfQueryResult SdfOracle::sparseContactAwareQuery(const Eigen::Vector3d& p) const
{
    return queryResultFromFast(queryFast(p, QueryMode::ContactAware));
}

// ---------- Witness estimation ----------

// ---------- Tricubic Hermite interpolation ----------
// Uses precomputed φ₀ and ∇φ₀ (n₀) at 2×2×2 = 8 corners.
// 1D Hermite basis:
//   h00(t) = 2t³ − 3t² + 1
//   h10(t) = t³ − 2t² + t
//   h01(t) = −2t³ + 3t²
//   h11(t) = t³ − t²
// Interpolant: p(t) = h00·f₀ + h10·f'₀ + h01·f₁ + h11·f'₁

SdfQueryResult SdfOracle::tricubicQuery(
    const Eigen::Vector3d& p) const
{
    const auto& spec = m_data.spec;

    // Fractional grid coordinates
    Eigen::Vector3d f = (p - spec.bmin).cwiseQuotient(spec.voxel)
        - Eigen::Vector3d(0.5, 0.5, 0.5);

    int i0 = static_cast<int>(std::floor(f.x()));
    int j0 = static_cast<int>(std::floor(f.y()));
    int k0 = static_cast<int>(std::floor(f.z()));

    double u = f.x() - i0;
    double v = f.y() - j0;
    double w = f.z() - k0;

    SdfQueryResult result;
    if (i0 < 0 || i0 + 1 >= spec.nx() ||
        j0 < 0 || j0 + 1 >= spec.ny() ||
        k0 < 0 || k0 + 1 >= spec.nz())
    {
        result.valid = false;
        return result;
    }

    // Gather φ and ∇φ at 8 corners
    double phi[2][2][2];
    double gx[2][2][2], gy[2][2][2], gz[2][2][2];
    double vx = spec.voxel.x(), vy = spec.voxel.y(), vz = spec.voxel.z();

    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        phi[di][dj][dk] = m_data.phi(i0 + di, j0 + dj, k0 + dk);
        const Eigen::Vector3d& n = m_data.normal(i0 + di, j0 + dj, k0 + dk);
        // ∇φ = n (unit gradient direction scaled by |∇φ|≈1 for SDF)
        gx[di][dj][dk] = n.x() * vx;  // ∂φ/∂u = ∂φ/∂x · Δx
        gy[di][dj][dk] = n.y() * vy;  // ∂φ/∂v = ∂φ/∂y · Δy
        gz[di][dj][dk] = n.z() * vz;  // ∂φ/∂w = ∂φ/∂z · Δz
    }

    // Hermite basis: h00, h10, h01, h11
    auto h00 = [](double t) { return t * t * (2.0 * t - 3.0) + 1.0; };
    auto h10 = [](double t) { return t * (t * (t - 2.0) + 1.0); };
    auto h01 = [](double t) { return t * t * (3.0 - 2.0 * t); };
    auto h11 = [](double t) { return t * t * (t - 1.0); };

    double a0 = h00(u), a1 = h10(u), a2 = h01(u), a3 = h11(u);
    double b0 = h00(v), b1 = h10(v), b2 = h01(v), b3 = h11(v);
    double c0 = h00(w), c1 = h10(w), c2 = h01(w), c3 = h11(w);

    // --- Interpolate φ via tensor-product Hermite ---
    double gap = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        double p0 = (di == 0) ? a0 : a2;
        double p1 = (di == 0) ? a1 : a3;
        double q0 = (dj == 0) ? b0 : b2;
        double q1 = (dj == 0) ? b1 : b3;
        double r0 = (dk == 0) ? c0 : c2;
        double r1 = (dk == 0) ? c1 : c3;

        // contribution: value + derivative terms
        gap += (p0 * q0 * r0) * phi[di][dj][dk];
        gap += (p1 * q0 * r0) * gx[di][dj][dk];
        gap += (p0 * q1 * r0) * gy[di][dj][dk];
        gap += (p0 * q0 * r1) * gz[di][dj][dk];
    }
    result.gap = gap;

    // --- Derivative basis functions ---
    // dh00/dt = 6t² − 6t,  dh10/dt = 3t² − 4t + 1
    // dh01/dt = −6t² + 6t, dh11/dt = 3t² − 2t
    auto dh00 = [](double t) { return 6.0 * t * (t - 1.0); };
    auto dh10 = [](double t) { return t * (3.0 * t - 4.0) + 1.0; };
    auto dh01 = [](double t) { return 6.0 * t * (1.0 - t); };
    auto dh11 = [](double t) { return t * (3.0 * t - 2.0); };

    double da0 = dh00(u), da1 = dh10(u), da2 = dh01(u), da3 = dh11(u);
    double db0 = dh00(v), db1 = dh10(v), db2 = dh01(v), db3 = dh11(v);
    double dc0 = dh00(w), dc1 = dh10(w), dc2 = dh01(w), dc3 = dh11(w);

    // ∂φ/∂u (then divide by vx for ∂φ/∂x)
    double dfdu = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        double p0 = (di == 0) ? a0 : a2;
        double p1 = (di == 0) ? a1 : a3;
        double dp0 = (di == 0) ? da0 : da2;
        double dp1 = (di == 0) ? da1 : da3;
        double q0 = (dj == 0) ? b0 : b2;
        double q1 = (dj == 0) ? b1 : b3;
        double r0 = (dk == 0) ? c0 : c2;
        double r1 = (dk == 0) ? c1 : c3;

        dfdu += (dp0 * q0 * r0) * phi[di][dj][dk];
        dfdu += (dp1 * q0 * r0) * gx[di][dj][dk];
        dfdu += (dp0 * q1 * r0) * gy[di][dj][dk];
        dfdu += (dp0 * q0 * r1) * gz[di][dj][dk];
    }

    // ∂φ/∂v
    double dfdv = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        double p0 = (di == 0) ? a0 : a2;
        double p1 = (di == 0) ? a1 : a3;
        double q0 = (dj == 0) ? b0 : b2;
        double q1 = (dj == 0) ? b1 : b3;
        double dq0 = (dj == 0) ? db0 : db2;
        double dq1 = (dj == 0) ? db1 : db3;
        double r0 = (dk == 0) ? c0 : c2;
        double r1 = (dk == 0) ? c1 : c3;

        dfdv += (p0 * dq0 * r0) * phi[di][dj][dk];
        dfdv += (p1 * dq0 * r0) * gx[di][dj][dk];
        dfdv += (p0 * dq1 * r0) * gy[di][dj][dk];
        dfdv += (p0 * dq0 * r1) * gz[di][dj][dk];
    }

    // ∂φ/∂w
    double dfdw = 0.0;
    for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
    for (int di = 0; di < 2; ++di)
    {
        double p0 = (di == 0) ? a0 : a2;
        double p1 = (di == 0) ? a1 : a3;
        double q0 = (dj == 0) ? b0 : b2;
        double q1 = (dj == 0) ? b1 : b3;
        double r0 = (dk == 0) ? c0 : c2;
        double r1 = (dk == 0) ? c1 : c3;
        double dr0 = (dk == 0) ? dc0 : dc2;
        double dr1 = (dk == 0) ? dc1 : dc3;

        dfdw += (p0 * q0 * dr0) * phi[di][dj][dk];
        dfdw += (p1 * q0 * dr0) * gx[di][dj][dk];
        dfdw += (p0 * q1 * dr0) * gy[di][dj][dk];
        dfdw += (p0 * q0 * dr1) * gz[di][dj][dk];
    }

    // Convert parametric derivatives to world-space gradient
    Eigen::Vector3d grad(dfdu / vx, dfdv / vy, dfdw / vz);
    double gLen = grad.norm();
    if (gLen > 1e-15)
        result.gapNormal = grad / gLen;
    else
        result.gapNormal = Eigen::Vector3d::UnitX();

    result.witnessPoint = p - result.gap * result.gapNormal;
    result.valid = true;
    return result;
}

// ---------- Witness estimation ----------

Eigen::Vector3d SdfOracle::estimateWitness(int i, int j, int k) const
{
    if (!m_data.spec.inBounds(i, j, k))
        return Eigen::Vector3d::Zero();
    return m_data.witness(i, j, k);
}

} // namespace sdf
