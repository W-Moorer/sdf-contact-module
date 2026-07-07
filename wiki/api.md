# Public API

## SdfOracle — Query Engine

```cpp
class SdfOracle {
    // Construction
    SdfOracle(const SdfData& data);
    SdfOracle(const SparseSdfData& data);
    static SdfOracle loadFromFile(const std::string& path);

    // Queries
    SdfQueryResult query(const Vector3d& point, QueryMode mode) const;
    SdfFastQueryResult queryFast(const Vector3d& point, QueryMode mode) const;
    std::vector<SdfQueryResult> queryBatch(
        const std::vector<Vector3d>& points, QueryMode mode) const;

    // Accessors
    const GridSpec& gridSpec() const;
    bool isValid() const;
    bool isSparse() const;
};
```

### Query Modes

| Mode | Description |
|------|-------------|
| `FirstOrder` | Taylor expansion: φ(p) = φ₀ + n₀·δx |
| `SecondOrder` | + Hessian: φ(p) = φ₀ + n₀·δx + ½δxᵀ·H₀·δx |
| `Trilinear` | Trilinear interpolation of φ₀ from 8 surrounding voxels |
| `Tricubic` | Hermite cubic interpolation using φ₀ and n₀ at 8 corners |
| `ContactAware` | Evaluates all precomputed contact candidates at the voxel |

### Query Results

```cpp
struct SdfQueryResult {
    double gap;                    // Signed distance
    Vector3d gapNormal;            // Unit outward normal
    Vector3d witnessPoint;        // Closest surface point
    bool valid;                   // Whether the query succeeded
    std::vector<SdfQueryCandidate> candidates;
};
```

## SdfIO — Serialization

```cpp
// Dense
void saveSdf(const SdfData& data, const std::string& path);
SdfData loadSdf(const std::string& path);

// Sparse
void saveSparseSdf(const SparseSdfData& data, const std::string& path);
SparseSdfData loadSparseSdf(const std::string& path);

// Validation
bool isValidSdf(const std::string& path);
bool isSparseSdf(const std::string& path);
```

## SdfBuilder — SDF Construction

```cpp
class SdfBuilder {
    static SdfData build(
        const std::vector<Vector3d>& vertices,
        const std::vector<Vector3i>& triangles,
        const BuilderParams& params);

    static SdfData build(
        const TriangleMesh& mesh,
        const BuilderParams& params);

    // Sparse
    static SparseSdfData buildSparse(
        const TriangleMesh& mesh,
        const SparseBuilderParams& params);

    // Contact-aware
    static SdfData buildContactAware(
        const CornerNormalTriangleMesh& mesh,
        const ContactAwareBuilderParams& params);

    static SparseSdfData buildSparseContactAware(
        const CornerNormalTriangleMesh& mesh,
        const ContactAwareBuilderParams& params);
};
```

## SdfSurfaceMesh — Mesh Loaders

```cpp
TriangleMesh loadObjFile(const std::string& path);
CornerNormalTriangleMesh loadContactAwareSurfaceFile(const std::string& path);
```

## NexDynAdapter

```cpp
using StandaloneSdfOracle = SdfOracle;
```
Thin alias for NexDyn integration.

Last verified against: `plugins/SdfOracle/include/SdfOracle/SdfOracle.h`, `plugins/SdfOracle/include/SdfOracle/SdfIO.h`, `plugins/SdfOracle/include/SdfOracle/SdfBuilder.h`, `plugins/SdfOracle/include/SdfOracle/SdfSurfaceMesh.h`
