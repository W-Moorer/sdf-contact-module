# plugins

C++ plugins for building, storing and querying Signed Distance Fields (SDF).

## Contents

| Directory | Description |
|-----------|-------------|
| `SdfOracle/` | Core SDF plugin: builder, query engine, I/O, demo |

## Dependencies

- Eigen3 (header-only, required)
- OpenMP (optional, for parallel acceleration)
- CUDA (optional, for GPU acceleration)

## Build

```bash
cmake -B build -S plugins/SdfOracle -DCMAKE_PREFIX_PATH="<path-to-eigen3>"
cmake --build build
```
