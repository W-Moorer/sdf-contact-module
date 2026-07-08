# Repository Layout

```
sdf-contact-module/
├── .opencode/
│   └── agents/
│       └── wiki-maintainer.md     # Agent for maintaining project wiki
│
├── assets/
│   └── validation/                # Comparison plots (gitignored)
│
├── docs/
│   ├── 多体动力学框架参考意见v1.md  # MBD framework reference v1 (Chinese)
│   └── 多体动力学框架参考意见v2.md  # MBD framework reference v2 — rank-deficient KKT (Chinese)
│
├── models/                        # 3D mesh & SDF files
│   ├── cube.obj / cube_res*.sdf   # Cube mesh & SDFs at various resolutions
│   ├── cylinder_res*.sdf          # Cylinder SDFs
│   ├── hollow_cylinder_*.obj      # Hollow cylinder meshes (multiple densities)
│   ├── rmd_*.obj                  # OBJs extracted from RecurDyn RMD GGEOM surfaces
│   └── (ReurDyn scenario dirs)    # RecurDyn model files (.rdyn, .rmd, .rplt)
│
├── plugins/
│   └── SdfOracle/
│       ├── CMakeLists.txt         # Build configuration
│       ├── include/               # Public headers (SdfCommon, GridSpec, SdfData, etc.)
│       ├── src/                   # Implementation (SdfOracle, SdfIO, SdfBuilder, etc.)
│       └── app/                   # SdfOracleDemo executable
│
├── prototype/
│   ├── src/
│   │   ├── ir/                    # MBD intermediate representation
│   │   │   ├── system_model.py    # SystemModel — top-level container
│   │   │   ├── body.py            # RigidBody
│   │   │   ├── frame.py           # Frame (marker/reference frame)
│   │   │   ├── joint.py           # Joint / FixedJoint / RevoluteJoint
│   │   │   ├── drive.py           # Drive / DriveMode
│   │   │   ├── contact.py         # ContactPair
│   │   │   ├── force.py           # ForceElement / TorqueElement
│   │   │   └── functions.py       # Function / Constant / Linear / Sin
│   │   │
│   │   ├── dynamics/              # MBD dynamics framework
│   │   │   ├── state.py           # State (r, R, v, ω)
│   │   │   ├── mass_matrix.py     # Block-diagonal mass matrix
│   │   │   ├── constraints.py     # Constraint residuals & Jacobians
│   │   │   ├── joints_fixed.py    # Fixed joint residual
│   │   │   ├── joints_revolute.py # Revolute joint residual & coordinate
│   │   │   ├── drives.py          # Motion drive residual
│   │   │   ├── contacts_sdf.py    # SDFContactEngine (vectorized, AABB broad+narrow phase)
│   │   │   ├── assembler.py       # DynamicsAssembler (M, J, Q, b_c)
│   │   │   └── utils.py           # Rotation utilities
│   │   │
│   │   ├── integrators/           # Time integrators
│   │   │   ├── explicit_projected.py  # Force-explicit + constraint KKT + projection
│   │   │   ├── implicit_be.py         # Backward Euler DAE (Newton + line search)
│   │   │   └── projection.py          # Pseudo-inverse position/velocity projection
│   │   │
│   │   ├── solvers/               # Linear solvers
│   │   │   └── kkt.py                 # Schur complement + SVD pseudo-inverse
│   │   │
│   │   ├── geometry/              # Contact geometry
│   │   │   ├── manifest.py        # GeometryManifest (body ↔ mesh/sdf mapping)
│   │   │   ├── quadrature.py      # QuadratureMesh, OBJ loading, centroid quadrature
│   │   │   └── sdf_cache.py       # SDF file cache
│   │   │
│   │   ├── importers/             # RMD file import
│   │   │   ├── rmd_lexer.py       # Tokenizer: KEYWORD / id → RMDEntity
│   │   │   ├── rmd_entities.py    # Raw data classes (RawBody, RawMarker, etc.)
│   │   │   ├── rmd_parser.py      # Parser: RMDEntity → RawRecurDynModel
│   │   │   ├── rmd_to_ir.py       # Normalizer: RawRecurDynModel → SystemModel
│   │   │   └── rmd_surface.py     # GGEOM surface extraction → OBJ via NODES/PATCHES
│   │   │
│   │   ├── validation/            # Validation framework
│   │   │   ├── signal_metrics.py  # RMSE/max error/smoke-strict classification
│   │   │   └── recur_compare.py   # RecurDynValidator, ValidationCase
│   │   │
│   │   ├── scenes/                # Scene definitions
│   │   ├── examples/              # 9 validation cases
│   │   ├── sdf_grid.py            # TrilinearSDFGrid, AnalyticPlaneSDF
│   │   ├── contact_law.py         # Normal pressure + friction traction
│   │   ├── contact_integrator.py  # Contact force integration
│   │   ├── mesh_quadrature.py     # Mesh loading & quadrature (legacy)
│   │   └── rigid_body.py          # RigidBody (legacy)
│   │
│   ├── data/                      # Sample RMD files
│   ├── run_level0.py              # Level 0–1 validation runner
│   └── run_mbd_examples.py        # MBD example runner
│
├── reference/
│   ├── rmd2obj.py                 # Reference: RMD GSurface → OBJ extraction
│   └── wiki-using.md              # Wiki maintenance instructions
│
├── scripts/                       # Utility scripts
│   ├── gen_meshes.py              # OBJ mesh generation (cube, hollow cylinder)
│   ├── export_all_signals.py      # Export RecurDyn RPLT → CSV via COM API
│   ├── compare_and_plot.py        # Three-way comparison (RecurDyn vs Framework vs Analytical)
│   └── export_from_recurdyn.py    # RecurDyn macro for in-application CSV export
│
├── wiki/                          # Project documentation
│   ├── index.md
│   ├── architecture.md
│   ├── api.md
│   ├── algorithms.md
│   ├── solver-backends.md
│   ├── validation.md
│   ├── build.md
│   └── layout.md
│
└── README.md
```

## Key Statistics

| Metric | Value |
|--------|-------|
| C++ source files | 6 (src/) |
| C++ header files | 10 (include/) |
| Total C++ lines | ~3200+ |
| Python modules | ~40+ files across 9 packages |
| Libraries | 2 static (SdfOracleCore, SdfOracle) |
| Executables | 1 (SdfOracleDemo) |
| Validation cases | 9 |
| Dependencies | Eigen3 (required), NumPy/SciPy (Python), RecurDyn 2023 (optional) |

Last verified against: repository file tree
