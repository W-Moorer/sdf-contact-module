// ============================================================================
// SdfBuilder.cpp — Build SdfData from triangle mesh
// ============================================================================

#include "SdfOracle/SdfBuilder.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sdf {

namespace {
using Clock = std::chrono::high_resolution_clock;

double elapsedMs(const Clock::time_point& start, const Clock::time_point& end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

Eigen::Vector3d normalizedOrDefault(
    const Eigen::Vector3d& value,
    const Eigen::Vector3d& fallback)
{
    const double norm = value.norm();
    if (norm <= 1.0e-15)
    {
        const double fallbackNorm = fallback.norm();
        if (fallbackNorm > 1.0e-15)
        {
            return fallback / fallbackNorm;
        }
        return Eigen::Vector3d::UnitX();
    }
    return value / norm;
}

Eigen::Vector3d barycentricOnTriangle(
    const Eigen::Vector3d& point,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c)
{
    const Eigen::Vector3d v0 = b - a;
    const Eigen::Vector3d v1 = c - a;
    const Eigen::Vector3d v2 = point - a;
    const double d00 = v0.dot(v0);
    const double d01 = v0.dot(v1);
    const double d11 = v1.dot(v1);
    const double d20 = v2.dot(v0);
    const double d21 = v2.dot(v1);
    const double denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) <= 1.0e-15)
    {
        return Eigen::Vector3d{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    }
    const double v = (d11 * d20 - d01 * d21) / denom;
    const double w = (d00 * d21 - d01 * d20) / denom;
    return Eigen::Vector3d{1.0 - v - w, v, w};
}

struct ContactAwareCellCandidate
{
    double phi{0.0};
    double distance{0.0};
    Eigen::Vector3d normal{1.0, 0.0, 0.0};
    Eigen::Vector3d witness{0.0, 0.0, 0.0};
    int branchId{-1};
};

struct ContactTriangleInfo
{
    Eigen::Vector3d a{0.0, 0.0, 0.0};
    Eigen::Vector3d b{0.0, 0.0, 0.0};
    Eigen::Vector3d c{0.0, 0.0, 0.0};
    std::array<Eigen::Vector3d, 3> cornerNormals{};
    Eigen::Vector3d faceNormal{1.0, 0.0, 0.0};
    Eigen::Vector3d bmin{0.0, 0.0, 0.0};
    Eigen::Vector3d bmax{0.0, 0.0, 0.0};
};

struct ContactTileIndex
{
    int tileDim{8};
    int nx{1};
    int ny{1};
    int nz{1};
    std::vector<std::vector<int>> triangleIds;

    int linearIndex(int ti, int tj, int tk) const
    {
        return ti + nx * (tj + ny * tk);
    }

    const std::vector<int>& candidatesForVoxel(int i, int j, int k) const
    {
        static const std::vector<int> empty;
        if (triangleIds.empty() || tileDim <= 0)
        {
            return empty;
        }
        const int ti = std::max(0, std::min(nx - 1, i / tileDim));
        const int tj = std::max(0, std::min(ny - 1, j / tileDim));
        const int tk = std::max(0, std::min(nz - 1, k / tileDim));
        return triangleIds[static_cast<size_t>(linearIndex(ti, tj, tk))];
    }
};

int branchIdFromNormal(const Eigen::Vector3d& rawNormal)
{
    const Eigen::Vector3d normal = normalizedOrDefault(rawNormal, Eigen::Vector3d::UnitX());
    auto quantize = [](double value) {
        value = std::max(-1.0, std::min(1.0, value));
        return static_cast<int32_t>(std::lround(value * 4095.0));
    };

    uint32_t hash = 2166136261u;
    const int32_t qx = quantize(normal.x());
    const int32_t qy = quantize(normal.y());
    const int32_t qz = quantize(normal.z());
    for (int32_t value : {qx, qy, qz})
    {
        hash ^= static_cast<uint32_t>(value);
        hash *= 16777619u;
    }
    return static_cast<int>(1u + (hash & 0x3fffffffu));
}

void mergeContactCandidate(
    std::vector<ContactAwareCellCandidate>& candidates,
    const ContactAwareCellCandidate& incoming,
    double cosMergeAngle)
{
    for (ContactAwareCellCandidate& existing : candidates)
    {
        if (existing.normal.dot(incoming.normal) >= cosMergeAngle)
        {
            if (std::abs(incoming.phi) < std::abs(existing.phi))
            {
                existing = incoming;
            }
            return;
        }
    }
    candidates.push_back(incoming);
}

std::vector<ContactTriangleInfo> makeContactTriangleInfos(
    const CornerNormalTriangleMesh& mesh)
{
    std::vector<ContactTriangleInfo> infos;
    infos.reserve(mesh.triangles.size());
    for (size_t triIndex = 0; triIndex < mesh.triangles.size(); ++triIndex)
    {
        const auto& tri = mesh.triangles[triIndex];
        for (int v : tri)
        {
            if (v < 0 || v >= static_cast<int>(mesh.vertices.size()))
            {
                throw std::runtime_error("Contact-aware mesh triangle references invalid vertex");
            }
        }

        ContactTriangleInfo info;
        info.a = mesh.vertices[tri[0]];
        info.b = mesh.vertices[tri[1]];
        info.c = mesh.vertices[tri[2]];
        info.cornerNormals = mesh.cornerNormals[triIndex];
        info.faceNormal = normalizedOrDefault(
            (info.b - info.a).cross(info.c - info.a),
            Eigen::Vector3d::UnitX());
        info.bmin = info.a.cwiseMin(info.b).cwiseMin(info.c);
        info.bmax = info.a.cwiseMax(info.b).cwiseMax(info.c);
        infos.push_back(info);
    }
    return infos;
}

int voxelCoordForWorldFloor(
    double value,
    double bmin,
    double voxel)
{
    return static_cast<int>(std::floor((value - bmin) / voxel - 0.5));
}

int voxelCoordForWorldCeil(
    double value,
    double bmin,
    double voxel)
{
    return static_cast<int>(std::ceil((value - bmin) / voxel - 0.5));
}

ContactTileIndex buildContactTileIndex(
    const GridSpec& spec,
    const std::vector<ContactTriangleInfo>& infos,
    int requestedTileDim,
    double influenceRadius)
{
    ContactTileIndex index;
    index.tileDim = std::max(1, requestedTileDim);
    index.nx = (spec.nx() + index.tileDim - 1) / index.tileDim;
    index.ny = (spec.ny() + index.tileDim - 1) / index.tileDim;
    index.nz = (spec.nz() + index.tileDim - 1) / index.tileDim;
    index.triangleIds.resize(static_cast<size_t>(index.nx * index.ny * index.nz));

    const double h = spec.voxel.maxCoeff();
    const double tileRadius =
        0.5 * std::sqrt(3.0) * static_cast<double>(index.tileDim) * h;
    const double expansion = std::max(0.0, influenceRadius) + tileRadius + h;

    for (int tri = 0; tri < static_cast<int>(infos.size()); ++tri)
    {
        const ContactTriangleInfo& info = infos[static_cast<size_t>(tri)];
        const Eigen::Vector3d lo =
            info.bmin - Eigen::Vector3d::Constant(expansion);
        const Eigen::Vector3d hi =
            info.bmax + Eigen::Vector3d::Constant(expansion);

        int i0 = voxelCoordForWorldFloor(lo.x(), spec.bmin.x(), spec.voxel.x());
        int j0 = voxelCoordForWorldFloor(lo.y(), spec.bmin.y(), spec.voxel.y());
        int k0 = voxelCoordForWorldFloor(lo.z(), spec.bmin.z(), spec.voxel.z());
        int i1 = voxelCoordForWorldCeil(hi.x(), spec.bmin.x(), spec.voxel.x());
        int j1 = voxelCoordForWorldCeil(hi.y(), spec.bmin.y(), spec.voxel.y());
        int k1 = voxelCoordForWorldCeil(hi.z(), spec.bmin.z(), spec.voxel.z());

        i0 = std::max(0, std::min(spec.nx() - 1, i0));
        j0 = std::max(0, std::min(spec.ny() - 1, j0));
        k0 = std::max(0, std::min(spec.nz() - 1, k0));
        i1 = std::max(0, std::min(spec.nx() - 1, i1));
        j1 = std::max(0, std::min(spec.ny() - 1, j1));
        k1 = std::max(0, std::min(spec.nz() - 1, k1));

        const int ti0 = i0 / index.tileDim;
        const int tj0 = j0 / index.tileDim;
        const int tk0 = k0 / index.tileDim;
        const int ti1 = i1 / index.tileDim;
        const int tj1 = j1 / index.tileDim;
        const int tk1 = k1 / index.tileDim;

        for (int tk = tk0; tk <= tk1; ++tk)
        for (int tj = tj0; tj <= tj1; ++tj)
        for (int ti = ti0; ti <= ti1; ++ti)
        {
            index.triangleIds[static_cast<size_t>(index.linearIndex(ti, tj, tk))]
                .push_back(tri);
        }
    }
    return index;
}

ContactAwareCellCandidate makeContactCandidate(
    const Eigen::Vector3d& point,
    const ContactTriangleInfo& info)
{
    const Eigen::Vector3d closest =
        closestPointOnTriangle(point, info.a, info.b, info.c);
    const Eigen::Vector3d bary =
        barycentricOnTriangle(closest, info.a, info.b, info.c);
    const Eigen::Vector3d normal = normalizedOrDefault(
        bary.x() * info.cornerNormals[0] +
            bary.y() * info.cornerNormals[1] +
            bary.z() * info.cornerNormals[2],
        info.faceNormal);
    return ContactAwareCellCandidate{
        (point - closest).dot(normal),
        (point - closest).norm(),
        normal,
        closest,
        branchIdFromNormal(normal)};
}

std::vector<ContactAwareCellCandidate> buildContactCandidatesForPoint(
    const Eigen::Vector3d& point,
    const std::vector<ContactTriangleInfo>& triangles,
    const std::vector<int>& candidateTriangleIds,
    const KdTree& nearestTree,
    double distanceTolerance,
    double cosMergeAngle,
    int maxCandidates)
{
    double minDistance = std::numeric_limits<double>::max();
    for (int triIndex : candidateTriangleIds)
    {
        if (triIndex < 0 || triIndex >= static_cast<int>(triangles.size()))
        {
            continue;
        }
        const ContactAwareCellCandidate candidate =
            makeContactCandidate(point, triangles[static_cast<size_t>(triIndex)]);
        minDistance = std::min(minDistance, candidate.distance);
    }

    std::array<int, 1> fallbackTriangleIds{-1};
    if (!std::isfinite(minDistance))
    {
        const ClosestPointResult nearest = nearestTree.findNearest(point);
        fallbackTriangleIds[0] = nearest.triangleIndex;
        minDistance = nearest.distance;
    }

    std::vector<ContactAwareCellCandidate> merged;
    merged.reserve(static_cast<size_t>(maxCandidates));
    const double allowedDistance = minDistance + distanceTolerance;
    auto tryTriangle = [&](int triIndex) {
        if (triIndex < 0 || triIndex >= static_cast<int>(triangles.size()))
        {
            return;
        }
        const ContactAwareCellCandidate candidate =
            makeContactCandidate(point, triangles[static_cast<size_t>(triIndex)]);
        if (candidate.distance <= allowedDistance)
        {
            mergeContactCandidate(merged, candidate, cosMergeAngle);
        }
    };

    if (candidateTriangleIds.empty())
    {
        tryTriangle(fallbackTriangleIds[0]);
    }
    else
    {
        for (int triIndex : candidateTriangleIds)
        {
            tryTriangle(triIndex);
        }
    }

    if (merged.empty())
    {
        if (fallbackTriangleIds[0] < 0)
        {
            const ClosestPointResult nearest = nearestTree.findNearest(point);
            fallbackTriangleIds[0] = nearest.triangleIndex;
        }
        tryTriangle(fallbackTriangleIds[0]);
    }

    std::sort(
        merged.begin(),
        merged.end(),
        [](const ContactAwareCellCandidate& a, const ContactAwareCellCandidate& b) {
            return std::abs(a.phi) < std::abs(b.phi);
        });
    if (static_cast<int>(merged.size()) > maxCandidates)
    {
        merged.resize(static_cast<size_t>(maxCandidates));
    }
    return merged;
}
} // namespace

// ---------- Public API ----------

SdfData SdfBuilder::build(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles,
    const BuilderParams& params)
{
    if (vertices.empty() || triangles.empty())
        throw std::runtime_error("Empty mesh provided to SdfBuilder");

    auto phaseStart = Clock::now();
    log(params, "[SdfBuilder] Setting up voxel grid...");
    GridSpec spec = setupGrid(vertices, params);
    log(params, "[SdfBuilder] Phase setup: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");
    log(params, "[SdfBuilder] Grid: " +
        std::to_string(spec.nx()) + " x " +
        std::to_string(spec.ny()) + " x " +
        std::to_string(spec.nz()) +
        " = " + std::to_string(spec.totalVoxels()) + " voxels");

    SdfData data;
    data.resize(spec);

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Building KD-tree...");
    KdTree kdTree;
    kdTree.build(vertices, triangles);
    log(params, "[SdfBuilder] Phase BVH build: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Computing SDF + normals...");
    computeSdfAndNormals(data, kdTree, vertices, triangles, params);
    log(params, "[SdfBuilder] Phase SDF fill: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Estimating Hessian...");
    estimateHessian(data);
    log(params, "[SdfBuilder] Phase Hessian: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    log(params, "[SdfBuilder] Done.");
    return data;
}

SdfData SdfBuilder::build(
    const TriangleMesh& mesh,
    const BuilderParams& params)
{
    return build(mesh.vertices, mesh.triangles, params);
}

SdfData SdfBuilder::buildContactAware(
    const CornerNormalTriangleMesh& mesh,
    const ContactAwareBuilderParams& params)
{
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::runtime_error("Empty mesh provided to contact-aware SdfBuilder");
    if (mesh.cornerNormals.size() != mesh.triangles.size())
        throw std::runtime_error("Contact-aware SDF requires one corner-normal triplet per triangle");
    if (params.maxCandidates <= 0)
        throw std::runtime_error("Contact-aware SDF maxCandidates must be positive");

    auto phaseStart = Clock::now();
    log(params, "[SdfBuilder] Setting up contact-aware voxel grid...");
    GridSpec spec = setupGrid(mesh.vertices, params);
    log(params, "[SdfBuilder] Phase setup: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");
    log(params, "[SdfBuilder] Grid: " +
        std::to_string(spec.nx()) + " x " +
        std::to_string(spec.ny()) + " x " +
        std::to_string(spec.nz()) +
        " = " + std::to_string(spec.totalVoxels()) + " voxels");

    SdfData data;
    data.resize(spec);
    data.resizeContactCandidates(params.maxCandidates);

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Preprocessing contact-aware triangles...");
    const std::vector<ContactTriangleInfo> contactTriangles =
        makeContactTriangleInfos(mesh);
    KdTree nearestTree;
    nearestTree.build(mesh.vertices, mesh.triangles);
    log(params, "[SdfBuilder] Phase contact-aware BVH build: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    const double h = spec.voxel.maxCoeff();
    const double distanceTolerance =
        params.featureDistanceTolerance > 0.0
            ? params.featureDistanceTolerance
            : 1.5 * h;
    const double mergeAngleRadians =
        params.normalMergeAngleDegrees * std::acos(-1.0) / 180.0;
    const double cosMergeAngle = std::cos(mergeAngleRadians);
    const int maxCandidates = data.contactCandidateStride;
    std::vector<int> allTriangleIds(contactTriangles.size());
    std::iota(allTriangleIds.begin(), allTriangleIds.end(), 0);
    const ContactTileIndex tileIndex =
        params.useTileAcceleration
            ? buildContactTileIndex(
                spec,
                contactTriangles,
                params.tileDim,
                distanceTolerance + 2.0 * h)
            : ContactTileIndex{};

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Computing contact-aware normal sectors...");
    const int total = spec.totalVoxels();
    const int progressStep = std::max(1, total / 20);
    int count = 0;

#pragma omp parallel for schedule(dynamic, 128)
    for (int idx = 0; idx < total; ++idx)
    {
        const int i = idx % spec.nx();
        const int j = (idx / spec.nx()) % spec.ny();
        const int k = idx / (spec.nx() * spec.ny());
        const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
        const std::vector<int>& tileCandidates =
            params.useTileAcceleration
                ? tileIndex.candidatesForVoxel(i, j, k)
                : allTriangleIds;
        const std::vector<ContactAwareCellCandidate> candidates =
            buildContactCandidatesForPoint(
                center,
                contactTriangles,
                tileCandidates,
                nearestTree,
                distanceTolerance,
                cosMergeAngle,
                maxCandidates);

        const int candidateCount =
            std::min(static_cast<int>(candidates.size()), maxCandidates);
        data.contactCandidateCount[static_cast<size_t>(idx)] =
            static_cast<uint8_t>(candidateCount);

        if (candidateCount > 0)
        {
            data.phi0[static_cast<size_t>(idx)] = candidates.front().phi;
            data.n0[static_cast<size_t>(idx)] = candidates.front().normal;
            data.witness0[static_cast<size_t>(idx)] = candidates.front().witness;
        }

        for (int candidate = 0; candidate < candidateCount; ++candidate)
        {
            const int out = data.contactCandidateIndex(idx, candidate);
            data.contactCandidatePhi0[static_cast<size_t>(out)] =
                candidates[static_cast<size_t>(candidate)].phi;
            data.contactCandidateNormal0[static_cast<size_t>(out)] =
                candidates[static_cast<size_t>(candidate)].normal;
            data.contactCandidateBranchId[static_cast<size_t>(out)] =
                candidates[static_cast<size_t>(candidate)].branchId;
        }

        if (params.logCallback)
        {
#pragma omp critical(SdfContactAwareBuilderProgress)
            {
                ++count;
                if (count % progressStep == 0)
                {
                    const int pct = (count * 100) / total;
                    params.logCallback("[SdfBuilder] Contact-aware progress: " +
                        std::to_string(pct) + "%");
                }
            }
        }
    }
    log(params, "[SdfBuilder] Phase contact-aware fill: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    if (params.storeHessian)
    {
        phaseStart = Clock::now();
        log(params, "[SdfBuilder] Estimating contact-aware primary Hessian...");
        estimateHessian(data);
        log(params, "[SdfBuilder] Phase Hessian: " +
            std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");
    }

    log(params, "[SdfBuilder] Done.");
    return data;
}

SparseSdfData SdfBuilder::buildSparseContactAware(
    const CornerNormalTriangleMesh& mesh,
    const ContactAwareBuilderParams& params)
{
    if (mesh.vertices.empty() || mesh.triangles.empty())
        throw std::runtime_error("Empty mesh provided to sparse contact-aware SdfBuilder");
    if (mesh.cornerNormals.size() != mesh.triangles.size())
        throw std::runtime_error("Sparse contact-aware SDF requires one corner-normal triplet per triangle");
    if (params.maxCandidates <= 0)
        throw std::runtime_error("Sparse contact-aware SDF maxCandidates must be positive");
    if (params.blockDim <= 0)
        throw std::runtime_error("Sparse contact-aware SDF blockDim must be positive");

    auto phaseStart = Clock::now();
    log(params, "[SdfBuilder] Setting up sparse contact-aware voxel grid...");
    GridSpec spec = setupGrid(mesh.vertices, params);
    const double h = spec.voxel.maxCoeff();
    const double bandDistance =
        params.bandDistance > 0.0
            ? params.bandDistance
            : static_cast<double>(params.bandVoxels) * h;
    const double activationDistance =
        bandDistance + static_cast<double>(params.haloVoxels) * h;

    SparseSdfData sparse;
    sparse.spec = spec;
    sparse.blockDim = params.blockDim;
    sparse.haloVoxels = params.haloVoxels;
    sparse.bandDistance = bandDistance;
    sparse.hasHessian = params.storeHessian;
    sparse.contactCandidateStride =
        std::max(0, std::min(params.maxCandidates, SdfData::kMaxContactCandidates));

    log(params, "[SdfBuilder] Grid: " +
        std::to_string(spec.nx()) + " x " +
        std::to_string(spec.ny()) + " x " +
        std::to_string(spec.nz()) +
        " = " + std::to_string(spec.totalVoxels()) + " dense voxels");
    log(params, "[SdfBuilder] Phase setup: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Preprocessing contact-aware triangles...");
    const std::vector<ContactTriangleInfo> contactTriangles =
        makeContactTriangleInfos(mesh);
    KdTree nearestTree;
    nearestTree.build(mesh.vertices, mesh.triangles);
    const double distanceTolerance =
        params.featureDistanceTolerance > 0.0
            ? params.featureDistanceTolerance
            : 1.5 * h;
    const ContactTileIndex tileIndex =
        params.useTileAcceleration
            ? buildContactTileIndex(
                spec,
                contactTriangles,
                params.tileDim,
                std::max(activationDistance, distanceTolerance + 2.0 * h))
            : ContactTileIndex{};
    log(params, "[SdfBuilder] Phase sparse contact-aware BVH/tile build: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    const double mergeAngleRadians =
        params.normalMergeAngleDegrees * std::acos(-1.0) / 180.0;
    const double cosMergeAngle = std::cos(mergeAngleRadians);
    const int maxCandidates = sparse.contactCandidateStride;
    std::vector<int> allTriangleIds(contactTriangles.size());
    std::iota(allTriangleIds.begin(), allTriangleIds.end(), 0);

    const int nbx = (spec.nx() + params.blockDim - 1) / params.blockDim;
    const int nby = (spec.ny() + params.blockDim - 1) / params.blockDim;
    const int nbz = (spec.nz() + params.blockDim - 1) / params.blockDim;
    const int totalBlocks = nbx * nby * nbz;
    const int blockVoxelCount = sparse.blockVoxelCount();

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Selecting sparse contact-aware blocks...");
    std::vector<std::uint8_t> active(static_cast<size_t>(totalBlocks), 0);

#pragma omp parallel for schedule(dynamic, 64)
    for (int blockLinear = 0; blockLinear < totalBlocks; ++blockLinear)
    {
        const int bi = blockLinear % nbx;
        const int bj = (blockLinear / nbx) % nby;
        const int bk = blockLinear / (nbx * nby);
        bool isActive = false;

        if (params.useTileAcceleration && !tileIndex.triangleIds.empty())
        {
            const int i0 = bi * params.blockDim;
            const int j0 = bj * params.blockDim;
            const int k0 = bk * params.blockDim;
            const int i1 = std::min(i0 + params.blockDim - 1, spec.nx() - 1);
            const int j1 = std::min(j0 + params.blockDim - 1, spec.ny() - 1);
            const int k1 = std::min(k0 + params.blockDim - 1, spec.nz() - 1);
            const int ti0 = i0 / tileIndex.tileDim;
            const int tj0 = j0 / tileIndex.tileDim;
            const int tk0 = k0 / tileIndex.tileDim;
            const int ti1 = i1 / tileIndex.tileDim;
            const int tj1 = j1 / tileIndex.tileDim;
            const int tk1 = k1 / tileIndex.tileDim;
            for (int tk = tk0; tk <= tk1 && !isActive; ++tk)
            for (int tj = tj0; tj <= tj1 && !isActive; ++tj)
            for (int ti = ti0; ti <= ti1; ++ti)
            {
                if (!tileIndex.triangleIds[
                        static_cast<size_t>(tileIndex.linearIndex(ti, tj, tk))]
                        .empty())
                {
                    isActive = true;
                    break;
                }
            }
        }
        else
        {
            const int i0 = bi * params.blockDim;
            const int j0 = bj * params.blockDim;
            const int k0 = bk * params.blockDim;
            const int i1 = std::min(i0 + params.blockDim, spec.nx());
            const int j1 = std::min(j0 + params.blockDim, spec.ny());
            const int k1 = std::min(k0 + params.blockDim, spec.nz());
            const Eigen::Vector3d blockCenter =
                spec.bmin + spec.voxel.cwiseProduct(Eigen::Vector3d{
                    0.5 * static_cast<double>(i0 + i1),
                    0.5 * static_cast<double>(j0 + j1),
                    0.5 * static_cast<double>(k0 + k1)});
            const ClosestPointResult nearest = nearestTree.findNearest(blockCenter);
            const double blockRadius =
                0.5 * std::sqrt(3.0) * static_cast<double>(params.blockDim) * h;
            isActive = nearest.distance <= activationDistance + blockRadius;
        }

        active[static_cast<size_t>(blockLinear)] = isActive ? 1 : 0;
    }

    for (int bk = 0; bk < nbz; ++bk)
    for (int bj = 0; bj < nby; ++bj)
    for (int bi = 0; bi < nbx; ++bi)
    {
        const int blockLinear = bi + nbx * (bj + nby * bk);
        if (active[static_cast<size_t>(blockLinear)] == 0)
        {
            continue;
        }

        SparseSdfBlock block;
        block.index = {bi, bj, bk};
        block.phi0.assign(blockVoxelCount, 0.0f);
        block.normal0.assign(3 * blockVoxelCount, 0.0f);
        if (params.storeHessian)
        {
            block.hessian0.assign(6 * blockVoxelCount, 0.0f);
        }
        block.contactCandidateCount.assign(blockVoxelCount, 0);
        block.contactCandidatePhi0.assign(
            static_cast<size_t>(blockVoxelCount * maxCandidates),
            0.0f);
        block.contactCandidateNormal0.assign(
            static_cast<size_t>(3 * blockVoxelCount * maxCandidates),
            0.0f);
        block.contactCandidateBranchId.assign(
            static_cast<size_t>(blockVoxelCount * maxCandidates),
            -1);
        sparse.blocks.emplace_back(std::move(block));
    }
    log(params, "[SdfBuilder] Phase sparse contact-aware block selection: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    std::sort(
        sparse.blocks.begin(),
        sparse.blocks.end(),
        [](const SparseSdfBlock& a, const SparseSdfBlock& b) {
            if (a.index[2] != b.index[2]) return a.index[2] < b.index[2];
            if (a.index[1] != b.index[1]) return a.index[1] < b.index[1];
            return a.index[0] < b.index[0];
        });

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Filling sparse contact-aware block samples...");
#pragma omp parallel for schedule(dynamic, 1)
    for (int blockIndex = 0; blockIndex < static_cast<int>(sparse.blocks.size()); ++blockIndex)
    {
        SparseSdfBlock& block = sparse.blocks[static_cast<size_t>(blockIndex)];
        for (int lz = 0; lz < params.blockDim; ++lz)
        for (int ly = 0; ly < params.blockDim; ++ly)
        for (int lx = 0; lx < params.blockDim; ++lx)
        {
            const int i = block.index[0] * params.blockDim + lx;
            const int j = block.index[1] * params.blockDim + ly;
            const int k = block.index[2] * params.blockDim + lz;
            const int local = block.localIndex(lx, ly, lz, params.blockDim);
            if (!spec.inBounds(i, j, k))
            {
                block.phi0[local] = 0.0f;
                block.normal0[3 * local + 0] = 1.0f;
                continue;
            }

            const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
            const std::vector<int>& tileCandidates =
                params.useTileAcceleration
                    ? tileIndex.candidatesForVoxel(i, j, k)
                    : allTriangleIds;
            const std::vector<ContactAwareCellCandidate> candidates =
                buildContactCandidatesForPoint(
                    center,
                    contactTriangles,
                    tileCandidates,
                    nearestTree,
                    distanceTolerance,
                    cosMergeAngle,
                    maxCandidates);
            const int candidateCount =
                std::min(static_cast<int>(candidates.size()), maxCandidates);
            block.contactCandidateCount[static_cast<size_t>(local)] =
                static_cast<uint8_t>(candidateCount);

            if (candidateCount > 0)
            {
                const ContactAwareCellCandidate& primary = candidates.front();
                block.phi0[local] = static_cast<float>(primary.phi);
                block.normal0[3 * local + 0] = static_cast<float>(primary.normal.x());
                block.normal0[3 * local + 1] = static_cast<float>(primary.normal.y());
                block.normal0[3 * local + 2] = static_cast<float>(primary.normal.z());
            }
            else
            {
                block.phi0[local] = 0.0f;
                block.normal0[3 * local + 0] = 1.0f;
                block.normal0[3 * local + 1] = 0.0f;
                block.normal0[3 * local + 2] = 0.0f;
            }

            for (int candidate = 0; candidate < candidateCount; ++candidate)
            {
                const int out =
                    block.contactCandidateIndex(local, candidate, maxCandidates);
                const ContactAwareCellCandidate& c =
                    candidates[static_cast<size_t>(candidate)];
                block.contactCandidatePhi0[static_cast<size_t>(out)] =
                    static_cast<float>(c.phi);
                block.contactCandidateNormal0[static_cast<size_t>(3 * out + 0)] =
                    static_cast<float>(c.normal.x());
                block.contactCandidateNormal0[static_cast<size_t>(3 * out + 1)] =
                    static_cast<float>(c.normal.y());
                block.contactCandidateNormal0[static_cast<size_t>(3 * out + 2)] =
                    static_cast<float>(c.normal.z());
                block.contactCandidateBranchId[static_cast<size_t>(out)] =
                    c.branchId;
            }
        }
    }
    log(params, "[SdfBuilder] Phase sparse contact-aware fill: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    sparse.rebuildIndex();
    if (params.storeHessian)
    {
        phaseStart = Clock::now();
        log(params, "[SdfBuilder] Estimating sparse contact-aware primary Hessian...");
        estimateSparseHessian(sparse);
        log(params, "[SdfBuilder] Phase sparse contact-aware Hessian: " +
            std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");
    }

    log(params, "[SdfBuilder] Sparse contact-aware SDF blocks: " +
        std::to_string(sparse.blocks.size()) + " / " +
        std::to_string(totalBlocks));
    return sparse;
}

SparseSdfData SdfBuilder::buildSparse(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles,
    const SparseBuilderParams& params)
{
    if (vertices.empty() || triangles.empty())
        throw std::runtime_error("Empty mesh provided to SdfBuilder");
    if (params.blockDim <= 0)
        throw std::runtime_error("Sparse SDF blockDim must be positive");

    log(params, "[SdfBuilder] Setting up sparse voxel grid...");
    GridSpec spec = setupGrid(vertices, params);
    const double h = spec.voxel.maxCoeff();
    const double bandDistance =
        params.bandDistance > 0.0 ? params.bandDistance
                                  : static_cast<double>(params.bandVoxels) * h;
    const double blockRadius =
        0.5 * std::sqrt(3.0) * static_cast<double>(params.blockDim) * h;
    const double activationDistance =
        bandDistance + static_cast<double>(params.haloVoxels) * h + blockRadius;

    auto phaseStart = Clock::now();
    SparseSdfData sparse;
    sparse.spec = spec;
    sparse.blockDim = params.blockDim;
    sparse.haloVoxels = params.haloVoxels;
    sparse.bandDistance = bandDistance;
    sparse.hasHessian = params.storeHessian;

    log(params, "[SdfBuilder] Building KD-tree...");
    KdTree kdTree;
    kdTree.build(vertices, triangles);
    log(params, "[SdfBuilder] Phase BVH build: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    const int nbx = (spec.nx() + params.blockDim - 1) / params.blockDim;
    const int nby = (spec.ny() + params.blockDim - 1) / params.blockDim;
    const int nbz = (spec.nz() + params.blockDim - 1) / params.blockDim;
    const int blockVoxelCount =
        params.blockDim * params.blockDim * params.blockDim;

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Selecting active sparse blocks...");
    const int totalBlocks = nbx * nby * nbz;
    std::vector<std::uint8_t> active(static_cast<size_t>(totalBlocks), 0);

#pragma omp parallel for schedule(dynamic, 64)
    for (int blockLinear = 0; blockLinear < totalBlocks; ++blockLinear)
    {
        const int bi = blockLinear % nbx;
        const int bj = (blockLinear / nbx) % nby;
        const int bk = blockLinear / (nbx * nby);
        const int i0 = bi * params.blockDim;
        const int j0 = bj * params.blockDim;
        const int k0 = bk * params.blockDim;
        const int i1 = std::min(i0 + params.blockDim, spec.nx());
        const int j1 = std::min(j0 + params.blockDim, spec.ny());
        const int k1 = std::min(k0 + params.blockDim, spec.nz());
        const Eigen::Vector3d blockCenter =
            spec.bmin + spec.voxel.cwiseProduct(Eigen::Vector3d{
                0.5 * static_cast<double>(i0 + i1),
                0.5 * static_cast<double>(j0 + j1),
                0.5 * static_cast<double>(k0 + k1)});

        const KdTree::SignedResult sr = kdTree.findSignedNearest(blockCenter);
        active[static_cast<size_t>(blockLinear)] =
            std::abs(sr.signedDistance) <= activationDistance ? 1 : 0;
    }

    for (int bk = 0; bk < nbz; ++bk)
    for (int bj = 0; bj < nby; ++bj)
    for (int bi = 0; bi < nbx; ++bi)
    {
        const int blockLinear = bi + nbx * (bj + nby * bk);
        if (active[static_cast<size_t>(blockLinear)] == 0)
        {
            continue;
        }

        SparseSdfBlock block;
        block.index = {bi, bj, bk};
        block.phi0.assign(blockVoxelCount, 0.0f);
        block.normal0.assign(3 * blockVoxelCount, 0.0f);
        if (params.storeHessian)
            block.hessian0.assign(6 * blockVoxelCount, 0.0f);
        sparse.blocks.emplace_back(std::move(block));
    }
    log(params, "[SdfBuilder] Phase sparse block selection: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    std::sort(
        sparse.blocks.begin(),
        sparse.blocks.end(),
        [](const SparseSdfBlock& a, const SparseSdfBlock& b) {
            if (a.index[2] != b.index[2]) return a.index[2] < b.index[2];
            if (a.index[1] != b.index[1]) return a.index[1] < b.index[1];
            return a.index[0] < b.index[0];
        });

    phaseStart = Clock::now();
    log(params, "[SdfBuilder] Filling sparse block samples...");
    const double eps = 1.0e-12;
#pragma omp parallel for schedule(dynamic, 1)
    for (int blockIndex = 0; blockIndex < static_cast<int>(sparse.blocks.size()); ++blockIndex)
    {
        SparseSdfBlock& block = sparse.blocks[static_cast<size_t>(blockIndex)];
        for (int lz = 0; lz < params.blockDim; ++lz)
        for (int ly = 0; ly < params.blockDim; ++ly)
        for (int lx = 0; lx < params.blockDim; ++lx)
        {
            const int i = block.index[0] * params.blockDim + lx;
            const int j = block.index[1] * params.blockDim + ly;
            const int k = block.index[2] * params.blockDim + lz;
            const int local = block.localIndex(lx, ly, lz, params.blockDim);
            if (!spec.inBounds(i, j, k))
            {
                block.phi0[local] = 0.0f;
                block.normal0[3 * local + 0] = 1.0f;
                continue;
            }

            const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);
            const KdTree::SignedResult sr = kdTree.findSignedNearest(center);
            const Eigen::Vector3d diff = center - sr.closestPoint;
            const double dist = diff.norm();
            Eigen::Vector3d grad;
            if (dist > eps)
            {
                const double sign = sr.signedDistance >= 0.0 ? 1.0 : -1.0;
                grad = sign * diff / dist;
            }
            else
            {
                const double sign = sr.signedDistance >= 0.0 ? 1.0 : -1.0;
                grad = sign * sr.normal;
            }

            block.phi0[local] = static_cast<float>(sr.signedDistance);
            block.normal0[3 * local + 0] = static_cast<float>(grad.x());
            block.normal0[3 * local + 1] = static_cast<float>(grad.y());
            block.normal0[3 * local + 2] = static_cast<float>(grad.z());
        }
    }
    log(params, "[SdfBuilder] Phase sparse fill: " +
        std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");

    sparse.rebuildIndex();
    if (params.storeHessian)
    {
        phaseStart = Clock::now();
        log(params, "[SdfBuilder] Estimating sparse Hessian...");
        estimateSparseHessian(sparse);
        log(params, "[SdfBuilder] Phase sparse Hessian: " +
            std::to_string(elapsedMs(phaseStart, Clock::now())) + " ms");
    }

    log(params, "[SdfBuilder] Sparse SDF blocks: " +
        std::to_string(sparse.blocks.size()) + " / " +
        std::to_string(nbx * nby * nbz));
    return sparse;
}

SparseSdfData SdfBuilder::buildSparse(
    const TriangleMesh& mesh,
    const SparseBuilderParams& params)
{
    return buildSparse(mesh.vertices, mesh.triangles, params);
}

// ---------- Grid setup ----------

GridSpec SdfBuilder::setupGrid(
    const std::vector<Eigen::Vector3d>& vertices,
    const BuilderParams& params)
{
    Eigen::Vector3d vmin(1e30, 1e30, 1e30);
    Eigen::Vector3d vmax(-1e30, -1e30, -1e30);
    for (const auto& v : vertices)
    {
        vmin = vmin.cwiseMin(v);
        vmax = vmax.cwiseMax(v);
    }

    double diag = (vmax - vmin).norm();
    double pad = params.padding * diag;
    Eigen::Vector3d bmin = vmin.array() - pad;
    Eigen::Vector3d bmax = vmax.array() + pad;

    Eigen::Vector3d size = bmax - bmin;
    double voxelSz = params.voxelSize;
    if (voxelSz <= 0.0)
    {
        double longest = size.maxCoeff();
        voxelSz = longest / static_cast<double>(params.targetResolution);
    }

    int rawNx = std::max(1, static_cast<int>(std::ceil(size.x() / voxelSz)));
    int rawNy = std::max(1, static_cast<int>(std::ceil(size.y() / voxelSz)));
    int rawNz = std::max(1, static_cast<int>(std::ceil(size.z() / voxelSz)));

    // Snap each dimension to standard power-of-2 resolution
    int nx = snapToStandardResolution(rawNx, params.maxResolution);
    int ny = snapToStandardResolution(rawNy, params.maxResolution);
    int nz = snapToStandardResolution(rawNz, params.maxResolution);

    GridSpec spec;
    spec.bmin = bmin;
    spec.voxel = Eigen::Vector3d(voxelSz, voxelSz, voxelSz);
    spec.shape = {nx, ny, nz};
    // Adjust bmax so that (bmax - bmin) / N = voxelSz exactly
    spec.bmax = bmin + Eigen::Vector3d(nx, ny, nz) * voxelSz;
    return spec;
}

// ---------- SDF + normal computation ----------

void SdfBuilder::computeSdfAndNormals(
    SdfData& data,
    const KdTree& kdTree,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::array<int, 3>>& triangles,
    const BuilderParams& params)
{
    const auto& spec = data.spec;
    const double eps = 1e-12;
    const int total = spec.totalVoxels();
    const int progressStep = std::max(1, total / 20);
    int count = 0;

#pragma omp parallel for schedule(dynamic, 256)
    for (int idx = 0; idx < total; ++idx)
    {
        const int i = idx % spec.nx();
        const int j = (idx / spec.nx()) % spec.ny();
        const int k = idx / (spec.nx() * spec.ny());
        const Eigen::Vector3d center = spec.centerFromIndex(i, j, k);

        // KD-tree nearest triangle query
        const KdTree::SignedResult sr = kdTree.findSignedNearest(center);

        // Store SDF value (outside positive, inside negative)
        data.phi0[static_cast<size_t>(idx)] = sr.signedDistance;

        // Store witness point
        data.witness0[static_cast<size_t>(idx)] = sr.closestPoint;

        // Compute gradient (unit outward normal)
        const Eigen::Vector3d diff = center - sr.closestPoint;
        const double dist = diff.norm();
        Eigen::Vector3d grad;
        if (dist > eps)
        {
            // sign(SDF) * (p - C) / ||p - C||
            const double sign = (sr.signedDistance >= 0.0) ? 1.0 : -1.0;
            grad = sign * diff / dist;
        }
        else
        {
            // Near surface: use face normal as fallback
            const double sign = (sr.signedDistance >= 0.0) ? 1.0 : -1.0;
            grad = sign * sr.normal;
        }
        data.n0[static_cast<size_t>(idx)] = grad;

        if (params.logCallback)
        {
#pragma omp critical(SdfBuilderProgress)
            {
                ++count;
                if (count % progressStep == 0)
                {
                    const int pct = (count * 100) / total;
                    params.logCallback("[SdfBuilder] SDF progress: " +
                        std::to_string(pct) + "%");
                }
            }
        }
    }
}

// ---------- Hessian estimation via central finite differences ----------

void SdfBuilder::estimateHessian(SdfData& data)
{
    const auto& spec = data.spec;
    const double vx = spec.voxel.x();
    const double vy = spec.voxel.y();
    const double vz = spec.voxel.z();

    // For each voxel, estimate Hessian columns from gradient differences.
    // H[:,axis] = (n[+1] - n[-1]) / (2 * spacing)
    // Then symmetrize: H = 0.5 * (H + H^T)

    // We need a temporary copy of normals to read neighbors safely.
    std::vector<Eigen::Vector3d> nCopy = data.n0;

    auto getNormal = [&](int i, int j, int k) -> Eigen::Vector3d
    {
        // Clamp to boundary
        i = std::max(0, std::min(i, spec.nx() - 1));
        j = std::max(0, std::min(j, spec.ny() - 1));
        k = std::max(0, std::min(k, spec.nz() - 1));
        return nCopy[spec.linearIndex(i, j, k)];
    };

    const int total = spec.totalVoxels();
#pragma omp parallel for schedule(static)
    for (int idx = 0; idx < total; ++idx)
    {
        const int i = idx % spec.nx();
        const int j = (idx / spec.nx()) % spec.ny();
        const int k = idx / (spec.nx() * spec.ny());
        // Central differences of normal field along each axis
        Eigen::Vector3d dn_dx = (getNormal(i+1, j, k) - getNormal(i-1, j, k)) / (2.0 * vx);
        Eigen::Vector3d dn_dy = (getNormal(i, j+1, k) - getNormal(i, j-1, k)) / (2.0 * vy);
        Eigen::Vector3d dn_dz = (getNormal(i, j, k+1) - getNormal(i, j, k-1)) / (2.0 * vz);

        // Build Hessian: H[row][col] = d(n[row])/d(x[col])
        // H = [dn_dx | dn_dy | dn_dz] as columns
        Eigen::Matrix3d H;
        H.col(0) = dn_dx;
        H.col(1) = dn_dy;
        H.col(2) = dn_dz;

        // Symmetrize
        H = 0.5 * (H + H.transpose());

        data.H0[static_cast<size_t>(idx)] = H;
    }
}

void SdfBuilder::estimateSparseHessian(SparseSdfData& data)
{
    const auto& spec = data.spec;
    const double vx = spec.voxel.x();
    const double vy = spec.voxel.y();
    const double vz = spec.voxel.z();

    auto normalOrFallback =
        [&](int i, int j, int k, const Eigen::Vector3d& fallback) {
            double phi = 0.0;
            Eigen::Vector3d normal = fallback;
            if (!data.sample(i, j, k, phi, normal, nullptr))
            {
                return fallback;
            }
            return normal;
        };

#pragma omp parallel for schedule(dynamic, 1)
    for (int blockIndex = 0; blockIndex < static_cast<int>(data.blocks.size()); ++blockIndex)
    {
        SparseSdfBlock& block = data.blocks[static_cast<size_t>(blockIndex)];
        if (block.hessian0.size() <
            static_cast<size_t>(6 * data.blockVoxelCount()))
        {
            block.hessian0.assign(6 * data.blockVoxelCount(), 0.0f);
        }

        for (int lz = 0; lz < data.blockDim; ++lz)
        for (int ly = 0; ly < data.blockDim; ++ly)
        for (int lx = 0; lx < data.blockDim; ++lx)
        {
            const int i = block.index[0] * data.blockDim + lx;
            const int j = block.index[1] * data.blockDim + ly;
            const int k = block.index[2] * data.blockDim + lz;
            if (!spec.inBounds(i, j, k))
            {
                continue;
            }

            const int local = block.localIndex(lx, ly, lz, data.blockDim);
            const Eigen::Vector3d current{
                static_cast<double>(block.normal0[3 * local + 0]),
                static_cast<double>(block.normal0[3 * local + 1]),
                static_cast<double>(block.normal0[3 * local + 2])};

            const Eigen::Vector3d dn_dx =
                (normalOrFallback(i + 1, j, k, current) -
                 normalOrFallback(i - 1, j, k, current)) / (2.0 * vx);
            const Eigen::Vector3d dn_dy =
                (normalOrFallback(i, j + 1, k, current) -
                 normalOrFallback(i, j - 1, k, current)) / (2.0 * vy);
            const Eigen::Vector3d dn_dz =
                (normalOrFallback(i, j, k + 1, current) -
                 normalOrFallback(i, j, k - 1, current)) / (2.0 * vz);

            Eigen::Matrix3d H;
            H.col(0) = dn_dx;
            H.col(1) = dn_dy;
            H.col(2) = dn_dz;
            H = 0.5 * (H + H.transpose());

            float* out = &block.hessian0[6 * local];
            out[0] = static_cast<float>(H(0, 0));
            out[1] = static_cast<float>(H(0, 1));
            out[2] = static_cast<float>(H(0, 2));
            out[3] = static_cast<float>(H(1, 1));
            out[4] = static_cast<float>(H(1, 2));
            out[5] = static_cast<float>(H(2, 2));
        }
    }
}

// ---------- Logging helper ----------

void SdfBuilder::log(const BuilderParams& params, const std::string& msg)
{
    if (params.logCallback)
        params.logCallback(msg);
    else
        std::cout << msg << std::endl;
}

} // namespace sdf
