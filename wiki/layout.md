# Repository Layout

```
sdf-contact-module/
├── .opencode/
│   └── agents/
│       └── wiki-maintainer.md     # Agent for maintaining project wiki
│
├── docs/
│   └── 参考意见v1.md               # Design notes — Python contact prototype (Chinese)
│
├── models/                        # 3D model files (OBJ, etc.) — currently empty
│
├── plugins/
│   └── SdfOracle/
│       ├── CMakeLists.txt         # Build configuration
│       ├── include/
│       │   └── SdfOracle/
│       │       ├── SdfCommon.h        # Query result types
│       │       ├── GridSpec.h         # Voxel grid specification
│       │       ├── SdfData.h          # Dense & sparse SDF data structures
│       │       ├── SdfOracle.h        # Query engine interface
│       │       ├── SdfIO.h            # Binary serialization interface
│       │       ├── SdfBuilder.h       # SDF builder interface
│       │       ├── SdfSurfaceMesh.h   # Mesh data structures
│       │       ├── KdTree.h           # BVH for nearest-triangle queries
│       │       └── NexDynAdapter.h    # NexDyn integration alias
│       ├── src/
│       │       ├── SdfCommon.cpp
│       │       ├── SdfOracle.cpp      # Query engine (~1100 lines)
│       │       ├── SdfIO.cpp          # Binary I/O (~554 lines)
│       │       ├── SdfBuilder.cpp     # Builder (~1237 lines)
│       │       └── SdfSurfaceMesh.cpp # Mesh loaders (~346 lines)
│       └── app/
│           └── SdfOracleDemo.cpp      # Standalone validation (~1007 lines)
│
├── reference/
│   └── wiki-using.md             # Wiki maintenance instructions
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
├── opencode.json                  # opencode project configuration
├── LICENSE
└── README.md
```

## Key Statistics

| Metric | Value |
|--------|-------|
| C++ source files | 6 (src/) |
| C++ header files | 10 (include/) |
| Total C++ lines | ~3200+ |
| Libraries | 2 static (SdfOracleCore, SdfOracle) |
| Executables | 1 (SdfOracleDemo) |
| Dependencies | Eigen3 (required), OpenMP (optional), CUDA (optional) |

Last verified against: repository file tree
