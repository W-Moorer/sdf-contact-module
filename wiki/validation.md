# Validation

## SdfOracleDemo

Standalone executable (`app/SdfOracleDemo.cpp`) for building, saving, loading, querying, and benchmarking SDFs.

### CLI Interface

```
SdfOracleDemo --obj <mesh.obj> --resolution <res> --out <prefix>
```

### Capabilities

| Feature | Description |
|---------|-------------|
| Build SDF | From OBJ with configurable resolution, padding, voxel size |
| Save/Load | Binary .sdf serialization round-trip |
| Query | Single-point and batch queries across all modes |
| Benchmark | Performance timing for query modes |
| CSV Export | Per-voxel diagnostic output |

## Python Prototype (Theory Validation)

`docs/参考意见v1.md` defines a Python prototype for contact integration validation with three levels:

### Level 0: Analytic Plane SDF
- Uses analytic plane `φ(x) = z - z_top` (no plugin SDF)
- Validates contact integration and friction formulas
- Expected: torsional friction torque matches analytic solution

### Level 1: Trilinear Cube SDF
- Plugin-generated cube `.sdf` at resolutions 64, 128, 256
- Python reads `.sdf` and performs trilinear queries
- Expected: results converge to Level 0 analytic solution as resolution increases

### Level 2: Full OBJ + SDF
- Uses actual `hollow_cylinder.obj` and `cube.obj`
- Validates mesh quadrature, SDF direction, gradient stability

### Acceptance Criteria

| Level | Metric | Threshold |
|-------|--------|-----------|
| 0 | Relative torque error | < 0.1% |
| 1 (res 64) | Relative torque error | < 5% |
| 1 (res 128) | Relative torque error | < 2% |
| 1 (res 256) | Relative torque error | < 1% |
| All | Tangential force residual | < 1-3% of μN |

Last verified against: `plugins/SdfOracle/app/SdfOracleDemo.cpp`, `docs/参考意见v1.md`
