# Build Workflows

## Dependencies

| Dependency | Version | Required | Notes |
|------------|---------|----------|-------|
| Eigen3 | 3.x | Yes | Header-only linear algebra |
| OpenMP | Any | No | Enabled by default, optional |
| CUDA | Any | No | Disabled by default, optional |
| CMake | ≥ 3.15 | Yes | |
| C++17 | ≥ 17 | Yes | |

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `NEXDYN_SDF_ENABLE_OPENMP` | `ON` | Enable OpenMP parallel acceleration |
| `NEXDYN_SDF_ENABLE_CUDA` | `OFF` | Enable CUDA GPU acceleration |

## Build Targets

| Target | Type | Description |
|--------|------|-------------|
| `SdfOracleCore` | Static library | Runtime query + I/O |
| `SdfOracle` | Static library | Core + Builder |
| `SdfOracleDemo` | Executable | Validation and benchmarking |

## Build Commands

```bash
# Configure (from project root)
cmake -B build -S plugins/SdfOracle

# With options
cmake -B build -S plugins/SdfOracle \
    -DNEXDYN_SDF_ENABLE_OPENMP=ON \
    -DNEXDYN_SDF_ENABLE_CUDA=OFF

# Build
cmake --build build

# Run demo
./build/SdfOracleDemo --obj model.obj --resolution 128 --out output
```

## MSVC-Specific

The CMake config sets these MSVC flags:
- `/utf-8` — UTF-8 source encoding
- `/wd4267 /wd4819 /wd4244` — Warning suppressions
- `/bigobj` — Extended object format
- `_CRT_SECURE_NO_WARNINGS` — Suppress secure CRT deprecation

Last verified against: `plugins/SdfOracle/CMakeLists.txt`
