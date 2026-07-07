#pragma once

// ============================================================================
// SdfBuilder — Build SdfData from triangle mesh
// ============================================================================
// Pipeline:
//   1. Setup voxel grid from mesh bounding box + padding/resolution
//   2. Build KD-tree from triangles
//   3. For each voxel: compute signed distance + normal + closest point
//   4. Estimate Hessian from normal field via central finite differences
//   5. Symmetrize Hessian: H = 0.5 * (H + H^T)
// ============================================================================

#include <functional>
#include <string>

#include "SdfOracle/SdfData.h"
#include "SdfOracle/KdTree.h"
#include "SdfOracle/SdfSurfaceMesh.h"

namespace sdf {

/// Parameters for SDF grid construction.
struct BuilderParams
{
    double padding{0.1};              ///< padding as fraction of bounding box diagonal
    double voxelSize{0.0};            ///< explicit voxel size (0 = auto from resolution)
    int targetResolution{128};        ///< target voxels along longest axis
    int maxResolution{512};           ///< maximum voxels along any axis
    std::function<void(const std::string&)> logCallback; ///< optional progress logger
};

struct SparseBuilderParams : public BuilderParams
{
    int blockDim{8};
    int haloVoxels{2};
    int bandVoxels{8};
    double bandDistance{0.0}; ///< explicit narrow-band half width (0 = bandVoxels * voxel)
    bool storeHessian{true};
};

struct ContactAwareBuilderParams : public BuilderParams
{
    int maxCandidates{4};
    double featureDistanceTolerance{0.0};
    double normalMergeAngleDegrees{5.0};
    int tileDim{8};
    bool useTileAcceleration{true};
    bool storeHessian{true};
    bool sparse{false};
    int blockDim{8};
    int haloVoxels{2};
    int bandVoxels{8};
    double bandDistance{0.0};
};

/// Build SdfData from a triangle mesh.
class SdfBuilder
{
public:
    /// Build from raw vertices + triangles.
    static SdfData build(
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<std::array<int, 3>>& triangles,
        const BuilderParams& params = {});

    /// Build from TriangleMesh struct.
    static SdfData build(
        const TriangleMesh& mesh,
        const BuilderParams& params = {});

    /// Build dense contact-aware SDF data from per-corner triangle normals.
    static SdfData buildContactAware(
        const CornerNormalTriangleMesh& mesh,
        const ContactAwareBuilderParams& params = {});

    /// Build a narrow-band sparse contact-aware SDF from per-corner triangle normals.
    static SparseSdfData buildSparseContactAware(
        const CornerNormalTriangleMesh& mesh,
        const ContactAwareBuilderParams& params = {});

    /// Build a narrow-band sparse SDF from raw vertices + triangles.
    static SparseSdfData buildSparse(
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<std::array<int, 3>>& triangles,
        const SparseBuilderParams& params = {});

    /// Build a narrow-band sparse SDF from TriangleMesh struct.
    static SparseSdfData buildSparse(
        const TriangleMesh& mesh,
        const SparseBuilderParams& params = {});

private:
    /// Setup grid spec from mesh bounds.
    static GridSpec setupGrid(
        const std::vector<Eigen::Vector3d>& vertices,
        const BuilderParams& params);

    /// Compute signed distance + normal + witness for all voxels.
    static void computeSdfAndNormals(
        SdfData& data,
        const KdTree& kdTree,
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<std::array<int, 3>>& triangles,
        const BuilderParams& params);

    /// Estimate Hessian from normal field via central finite differences.
    static void estimateHessian(SdfData& data);

    static void estimateSparseHessian(SparseSdfData& data);

    static void log(const BuilderParams& params, const std::string& msg);
};

} // namespace sdf
