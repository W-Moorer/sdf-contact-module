#pragma once

// ============================================================================
// SdfData — Dense voxel data container
// ============================================================================
// Stores per-voxel:
//   phi0  : SDF value at voxel center
//   n0    : unit outward normal (∇SDF) at voxel center
//   H0    : 3×3 symmetric Hessian matrix (6 independent components)
//
// Storage layout: flat arrays indexed by linearIndex(i,j,k).
// Hessian stored as 6 components: Hxx, Hxy, Hxz, Hyy, Hyz, Hzz
// ============================================================================

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "SdfOracle/GridSpec.h"
#include "SdfOracle/SdfCommon.h"

namespace sdf {

struct SdfData
{
    static constexpr int kMaxContactCandidates = kMaxSdfQueryCandidates;

    GridSpec spec;

    /// Per-voxel SDF value, shape = totalVoxels
    std::vector<double> phi0;

    /// Per-voxel unit outward normal, shape = totalVoxels × 3
    std::vector<Eigen::Vector3d> n0;

    /// Per-voxel Hessian (symmetric 3×3), stored as 6 components
    /// [Hxx, Hxy, Hxz, Hyy, Hyz, Hzz], shape = totalVoxels
    std::vector<Eigen::Matrix3d> H0;

    /// Per-voxel closest surface point (for witness point estimation)
    std::vector<Eigen::Vector3d> witness0;

    int contactCandidateStride{0};
    std::vector<uint8_t> contactCandidateCount;
    std::vector<double> contactCandidatePhi0;
    std::vector<Eigen::Vector3d> contactCandidateNormal0;
    std::vector<int32_t> contactCandidateBranchId;

    void resize(const GridSpec& s)
    {
        spec = s;
        int n = s.totalVoxels();
        phi0.assign(n, 0.0);
        n0.assign(n, Eigen::Vector3d::UnitX());
        H0.assign(n, Eigen::Matrix3d::Zero());
        witness0.assign(n, Eigen::Vector3d::Zero());
        clearContactCandidates();
    }

    void resizeContactCandidates(int maxCandidates)
    {
        const int stride = std::max(0, std::min(maxCandidates, kMaxContactCandidates));
        contactCandidateStride = stride;
        const int n = spec.totalVoxels();
        contactCandidateCount.assign(n, 0);
        contactCandidatePhi0.assign(static_cast<size_t>(n * stride), 0.0);
        contactCandidateNormal0.assign(
            static_cast<size_t>(n * stride),
            Eigen::Vector3d::UnitX());
        contactCandidateBranchId.assign(static_cast<size_t>(n * stride), -1);
    }

    void clearContactCandidates()
    {
        contactCandidateStride = 0;
        contactCandidateCount.clear();
        contactCandidatePhi0.clear();
        contactCandidateNormal0.clear();
        contactCandidateBranchId.clear();
    }

    bool hasContactCandidates() const
    {
        const int n = spec.totalVoxels();
        return contactCandidateStride > 0 &&
            contactCandidateCount.size() == static_cast<size_t>(n) &&
            contactCandidatePhi0.size() ==
                static_cast<size_t>(n * contactCandidateStride) &&
            contactCandidateNormal0.size() ==
                static_cast<size_t>(n * contactCandidateStride) &&
            contactCandidateBranchId.size() ==
                static_cast<size_t>(n * contactCandidateStride);
    }

    int contactCandidateIndex(int voxelIndex, int candidateIndex) const
    {
        return voxelIndex * contactCandidateStride + candidateIndex;
    }

    /// Access phi0 at (i,j,k)
    double& phi(int i, int j, int k) { return phi0[spec.linearIndex(i, j, k)]; }
    double phi(int i, int j, int k) const { return phi0[spec.linearIndex(i, j, k)]; }

    /// Access n0 at (i,j,k)
    Eigen::Vector3d& normal(int i, int j, int k) { return n0[spec.linearIndex(i, j, k)]; }
    const Eigen::Vector3d& normal(int i, int j, int k) const { return n0[spec.linearIndex(i, j, k)]; }

    /// Access H0 at (i,j,k)
    Eigen::Matrix3d& hessian(int i, int j, int k) { return H0[spec.linearIndex(i, j, k)]; }
    const Eigen::Matrix3d& hessian(int i, int j, int k) const { return H0[spec.linearIndex(i, j, k)]; }

    /// Access witness at (i,j,k)
    Eigen::Vector3d& witness(int i, int j, int k) { return witness0[spec.linearIndex(i, j, k)]; }
    const Eigen::Vector3d& witness(int i, int j, int k) const { return witness0[spec.linearIndex(i, j, k)]; }
};

struct SparseSdfBlock
{
    std::array<int, 3> index{0, 0, 0};
    std::vector<float> phi0;
    std::vector<float> normal0;  // x,y,z per local voxel
    std::vector<float> hessian0; // Hxx,Hxy,Hxz,Hyy,Hyz,Hzz per local voxel
    std::vector<uint8_t> contactCandidateCount;
    std::vector<float> contactCandidatePhi0;
    std::vector<float> contactCandidateNormal0; // x,y,z per local candidate
    std::vector<int32_t> contactCandidateBranchId;

    int localIndex(int x, int y, int z, int blockDim) const
    {
        return x + blockDim * (y + blockDim * z);
    }

    int contactCandidateIndex(int localVoxelIndex, int candidateIndex, int stride) const
    {
        return localVoxelIndex * stride + candidateIndex;
    }

    bool hasContactCandidates(int blockVoxelCount, int stride) const
    {
        return stride > 0 &&
            contactCandidateCount.size() >= static_cast<size_t>(blockVoxelCount) &&
            contactCandidatePhi0.size() >=
                static_cast<size_t>(blockVoxelCount * stride) &&
            contactCandidateNormal0.size() >=
                static_cast<size_t>(3 * blockVoxelCount * stride) &&
            contactCandidateBranchId.size() >=
                static_cast<size_t>(blockVoxelCount * stride);
    }
};

struct SparseSdfData
{
    GridSpec spec;
    int blockDim{8};
    int haloVoxels{2};
    double bandDistance{0.0};
    bool hasHessian{true};
    int contactCandidateStride{0};
    std::vector<SparseSdfBlock> blocks;
    std::unordered_map<uint64_t, size_t> blockLookup;

    static uint64_t key(int bi, int bj, int bk)
    {
        return static_cast<uint64_t>(static_cast<uint32_t>(bi)) |
            (static_cast<uint64_t>(static_cast<uint32_t>(bj)) << 21) |
            (static_cast<uint64_t>(static_cast<uint32_t>(bk)) << 42);
    }

    void rebuildIndex()
    {
        blockLookup.clear();
        blockLookup.reserve(blocks.size());
        for (size_t i = 0; i < blocks.size(); ++i)
        {
            blockLookup.emplace(
                key(blocks[i].index[0], blocks[i].index[1], blocks[i].index[2]),
                i);
        }
    }

    bool empty() const
    {
        return blocks.empty();
    }

    int blockVoxelCount() const
    {
        return blockDim * blockDim * blockDim;
    }

    size_t totalStoredVoxels() const
    {
        return blocks.size() * static_cast<size_t>(blockVoxelCount());
    }

    bool hasContactCandidates() const
    {
        if (contactCandidateStride <= 0 || blocks.empty())
        {
            return false;
        }
        const int voxelCount = blockVoxelCount();
        for (const SparseSdfBlock& block : blocks)
        {
            if (!block.hasContactCandidates(voxelCount, contactCandidateStride))
            {
                return false;
            }
        }
        return true;
    }

    const SparseSdfBlock* blockForVoxel(int i, int j, int k) const
    {
        if (!spec.inBounds(i, j, k) || blockDim <= 0)
        {
            return nullptr;
        }
        const int bi = i / blockDim;
        const int bj = j / blockDim;
        const int bk = k / blockDim;
        const auto it = blockLookup.find(key(bi, bj, bk));
        if (it == blockLookup.end())
        {
            return nullptr;
        }
        return &blocks[it->second];
    }

    bool sample(
        int i,
        int j,
        int k,
        double& phi,
        Eigen::Vector3d& normal,
        Eigen::Matrix3d* hessian = nullptr) const
    {
        const SparseSdfBlock* block = blockForVoxel(i, j, k);
        if (block == nullptr)
        {
            return false;
        }
        return sampleFromBlock(*block, i, j, k, phi, normal, hessian);
    }

    bool sampleFromBlock(
        const SparseSdfBlock& block,
        int i,
        int j,
        int k,
        double& phi,
        Eigen::Vector3d& normal,
        Eigen::Matrix3d* hessian = nullptr) const
    {
        const int lx = i - block.index[0] * blockDim;
        const int ly = j - block.index[1] * blockDim;
        const int lz = k - block.index[2] * blockDim;
        if (lx < 0 || lx >= blockDim || ly < 0 || ly >= blockDim ||
            lz < 0 || lz >= blockDim)
        {
            return false;
        }

        const int local = block.localIndex(lx, ly, lz, blockDim);
        const int voxelCount = blockVoxelCount();
        if (local < 0 || local >= voxelCount ||
            block.phi0.size() < static_cast<size_t>(voxelCount) ||
            block.normal0.size() < static_cast<size_t>(3 * voxelCount))
        {
            return false;
        }

        phi = static_cast<double>(block.phi0[local]);
        normal = Eigen::Vector3d{
            static_cast<double>(block.normal0[3 * local + 0]),
            static_cast<double>(block.normal0[3 * local + 1]),
            static_cast<double>(block.normal0[3 * local + 2])};

        if (hessian != nullptr)
        {
            if (!hasHessian ||
                block.hessian0.size() < static_cast<size_t>(6 * voxelCount))
            {
                return false;
            }
            const float* h = &block.hessian0[6 * local];
            (*hessian)(0, 0) = h[0];
            (*hessian)(0, 1) = h[1];
            (*hessian)(0, 2) = h[2];
            (*hessian)(1, 0) = h[1];
            (*hessian)(1, 1) = h[3];
            (*hessian)(1, 2) = h[4];
            (*hessian)(2, 0) = h[2];
            (*hessian)(2, 1) = h[4];
            (*hessian)(2, 2) = h[5];
        }

        return true;
    }
};

} // namespace sdf
