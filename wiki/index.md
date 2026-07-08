# sdf-contact-module Wiki

A **RecurDyn RMD-compatible multibody dynamics validation framework** with **SDF (Signed Distance Field) contact mechanics**.

## Project Components

| Layer | Description | Language |
|-------|-------------|----------|
| **SdfOracle** | Build, store & query SDF from triangle meshes | C++17 |
| **MBD-IR** | Multibody dynamics intermediate representation | Python |
| **Dynamics** | Maximal-coordinate constrained dynamics (KKT + SVD) | Python |
| **Contact** | AABB-BVH broad phase + trilinear SDF narrow phase | Python |
| **Importers** | RecurDyn RMD format parser (real format) | Python |
| **Validation** | RecurDyn vs Framework vs Analytical comparison | Python |

## Wiki Pages

| Page | Description |
|------|-------------|
| [Architecture](architecture.md) | SdfOracle plugin layering (C++) |
| [Public API](api.md) | Query engine, I/O, builder interfaces |
| [Algorithms](algorithms.md) | SDF query modes: trilinear, tricubic, contact-aware |
| [Solver Backends](solver-backends.md) | OpenMP and CUDA acceleration |
| [Validation](validation.md) | Validation cases & comparison plots |
| [Build](build.md) | CMake configuration and dependencies |
| [Repository Layout](layout.md) | Directory structure and file organization |

## Quick Links

- **Python MBD framework**: `prototype/src/`
- **SDF C++ plugin**: `plugins/SdfOracle/`
- **RecurDyn RMD importers**: `prototype/src/importers/`
- **Validation scripts**: `scripts/compare_and_plot.py`
- **Model files**: `models/` (OBJ + SDF + RMD)

Last verified against: repository file tree
