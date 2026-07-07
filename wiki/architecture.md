# Architecture

SdfOracle is organized in three layers, each building on the one below.

## Layer 1: Core (`SdfOracleCore`)

Static library with runtime-only dependencies. No mesh or builder code.

| Component | File | Responsibility |
|-----------|------|----------------|
| `GridSpec` | `GridSpec.h` | Cell-centered voxel grid definition; coordinate transforms (world ↔ voxel index); resolution snapping |
| `SdfData` / `SparseSdfData` | `SdfData.h` | Dense / sparse narrow-band SDF storage containers with per-voxel phi, normal, Hessian, witness, and contact candidates |
| `SdfOracle` | `SdfOracle.h/.cpp` | Query engine dispatching to five modes for dense and sparse data |
| `SdfIO` | `SdfIO.h/.cpp` | Binary serialization of SDF data to/from `.sdf` files |

## Layer 2: Builder (`SdfOracle`)

Static library extending core with mesh processing capabilities.

| Component | File | Responsibility |
|-----------|------|----------------|
| `SdfBuilder` | `SdfBuilder.h/.cpp` | Builds dense/sparse/contact-aware SDF from triangle meshes; Hessian estimation |
| `SdfSurfaceMesh` | `SdfSurfaceMesh.h/.cpp` | Mesh loaders for OBJ and custom contact-aware surface format |
| `KdTree` | `KdTree.h` | Header-only BVH for nearest-triangle queries and inside/outside test |

## Layer 3: Demo (`SdfOracleDemo`)

Standalone executable for validation.

| Component | File | Responsibility |
|-----------|------|----------------|
| `SdfOracleDemo` | `app/SdfOracleDemo.cpp` | CLI interface: build, save, load, query, benchmark, export CSV |

## Data Flow

```
OBJ file
   → SdfSurfaceMesh::loadObjFile()
   → TriangleMesh (vertices + triangles)
   → SdfBuilder::build(mesh, params)
     → KdTree::build(vertices, triangles)
     → For each voxel center: KdTree::findSignedNearest() → phi0, n0, witness0
     → SdfBuilder::estimateHessian() → H0 (central differences on normal field)
   → SdfData (phi0, n0, H0, witness0)
     → SdfIO::saveSdf() → .sdf binary file
     → SdfOracle(data) → query engine
       → query(point, mode) → SdfQueryResult { gap, gapNormal, witnessPoint }
```

## Design Decisions

- **Cell-centered grid**: Values stored at voxel centers with +0.5 offset in index transforms
- **Dual precision**: Dense uses `double`; sparse uses `float` for memory efficiency
- **Power-of-2 resolution snapping**: Grid dimensions snap to standard resolutions `{8, 16, 32, 64, 128, 256, 512, 1024, 2048}`
- **Sparse narrow-band**: Only voxels within a band around the surface are stored, using 8×8×8 blocks
- **Contact-aware extension**: Per-voxel multiple surface-feature candidates for multi-contact resolution

Last verified against: `plugins/SdfOracle/CMakeLists.txt`, `plugins/SdfOracle/include/SdfOracle/*.h`, `plugins/SdfOracle/src/*.cpp`
