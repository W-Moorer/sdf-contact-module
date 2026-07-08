# prototype

Python multibody dynamics validation framework with SDF contact mechanics.

## Quick Start

```bash
cd prototype
python run_mbd_examples.py    # All validation cases
python run_level0.py          # Original contact benchmark
```

## Layout

```
src/
├── ir/           MBD-IR intermediate representation
├── dynamics/     Multibody dynamics (constraints, mass matrix, KKT)
├── integrators/  Explicit projected + implicit Backward Euler
├── solvers/      KKT solver with SVD pseudo-inverse
├── geometry/     Quadrature mesh, SDF cache, geometry manifest
├── importers/    RecurDyn RMD format parser
├── validation/   Validation framework (metrics, comparison)
├── scenes/       Scene definitions
└── examples/     9 validation cases
```
