#pragma once

// ============================================================================
// GridSpec — Voxel grid specification and coordinate transforms
// ============================================================================
// Cell-centered grid: voxel center = bmin + voxel * (ijk + 0.5)
// World → index:    ijk = round((p - bmin) / voxel - 0.5)
// ============================================================================

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <tuple>

namespace sdf {

struct GridSpec
{
    Eigen::Vector3d bmin{0.0, 0.0, 0.0};
    Eigen::Vector3d bmax{1.0, 1.0, 1.0};
    Eigen::Vector3d voxel{1.0, 1.0, 1.0};
    std::array<int, 3> shape{1, 1, 1}; // (nx, ny, nz)

    int nx() const { return shape[0]; }
    int ny() const { return shape[1]; }
    int nz() const { return shape[2]; }
    int totalVoxels() const { return shape[0] * shape[1] * shape[2]; }

    /// Voxel center from grid index.
    Eigen::Vector3d centerFromIndex(int i, int j, int k) const
    {
        return bmin + voxel.cwiseProduct(
            Eigen::Vector3d(i + 0.5, j + 0.5, k + 0.5));
    }

    /// Nearest voxel index from world coordinate (round-to-nearest).
    void indexFromWorld(const Eigen::Vector3d& p,
                        int& i, int& j, int& k) const
    {
        Eigen::Vector3d f = (p - bmin).cwiseQuotient(voxel)
            - Eigen::Vector3d(0.5, 0.5, 0.5);
        i = static_cast<int>(std::round(f.x()));
        j = static_cast<int>(std::round(f.y()));
        k = static_cast<int>(std::round(f.z()));
    }

    /// Check if voxel index is within bounds.
    bool inBounds(int i, int j, int k) const
    {
        return i >= 0 && i < shape[0]
            && j >= 0 && j < shape[1]
            && k >= 0 && k < shape[2];
    }

    /// Linear index from 3D index.
    int linearIndex(int i, int j, int k) const
    {
        return i + shape[0] * (j + shape[1] * k);
    }
};

/// Snap dimension to nearest standard resolution (power-of-2 ladder).
/// Finds the smallest value in {8,16,32,64,128,256,512,1024,2048} that is >= n,
/// clamped to maxRes.
inline int snapToStandardResolution(int n, int maxRes)
{
    static constexpr int levels[] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    for (int v : levels)
    {
        if (v >= n || v >= maxRes)
            return std::min(v, maxRes);
    }
    return maxRes;
}

} // namespace sdf
