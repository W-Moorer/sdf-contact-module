# sdf-contact-module Wiki

SdfOracle — a standalone C++ plugin for building, storing, and querying Signed Distance Fields (SDF) from triangle meshes. Designed for contact mechanics in rigid-body simulation.

## Wiki Pages

| Page | Description |
|------|-------------|
| [Architecture](architecture.md) | Overall plugin architecture and layering |
| [Public API](api.md) | Query engine, I/O, builder interfaces |
| [Algorithms](algorithms.md) | Query modes: trilinear, first/second-order, tricubic, contact-aware |
| [Solver Backends](solver-backends.md) | OpenMP and CUDA acceleration |
| [Validation](validation.md) | SdfOracleDemo executable, CLI, benchmarks |
| [Build](build.md) | CMake configuration and dependencies |
| [Repository Layout](layout.md) | Directory structure and file organization |

## Quick Links

- **Source**: `plugins/SdfOracle/`
- **Build target**: `SdfOracleCore` (static lib), `SdfOracle` (static lib), `SdfOracleDemo` (executable)
- **Dependencies**: Eigen3 (required), OpenMP (optional), CUDA (optional)
- **File format**: `.sdf` binary (magic `"SDFO"`, version 0.1.0)

Last verified against: `plugins/SdfOracle/CMakeLists.txt`, `plugins/SdfOracle/include/SdfOracle/*.h`
