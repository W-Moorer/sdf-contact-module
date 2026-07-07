# Solver Backends

## OpenMP (Default: ON)

Controlled by `NEXDYN_SDF_ENABLE_OPENMP` CMake option.

- Used in `SdfBuilder::computeSdfAndNormals()` — parallel over voxels
- Used in `SdfBuilder::estimateHessian()` — parallel over voxels
- Used in `SdfOracle::queryBatch()` — parallel over query points
- Guarded by `SDF_ORACLE_HAS_OPENMP` compile definition

```cmake
option(NEXDYN_SDF_ENABLE_OPENMP "Enable OpenMP acceleration" ON)
find_package(OpenMP QUIET)
```

## CUDA (Default: OFF)

Controlled by `NEXDYN_SDF_ENABLE_CUDA` CMake option.

- Enabled only when CUDA compiler is available
- Guarded by `SDF_ORACLE_ENABLE_CUDA` compile definition
- Links `CUDA::cudart`

```cmake
option(NEXDYN_SDF_ENABLE_CUDA "Enable optional CUDA acceleration" OFF)
check_language(CUDA)
```

## Eigen3 (Required)

Header-only linear algebra library. Used throughout for `Vector3d`, `Matrix3d`, and related operations. No version constraint beyond C++17 compatibility.

Last verified against: `plugins/SdfOracle/CMakeLists.txt`
