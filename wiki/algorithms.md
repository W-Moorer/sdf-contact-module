# Contact Algorithms

## Query Modes

### First Order

Taylor expansion at nearest voxel center:

```
φ(p) = φ₀ + n₀ · (p - c)
```

where `c` is the nearest voxel center, `φ₀` is the stored SDF value, and `n₀` is the stored unit normal.

### Second Order

Adds Hessian term for improved accuracy away from voxel centers:

```
φ(p) = φ₀ + n₀ · δx + 0.5 · δxᵀ · H₀ · δx
```

Hessian `H₀` is estimated via central finite differences on the normal field, then symmetrized.

### Trilinear (Primary Contact Mode)

Trilinear interpolation of `φ₀` from the 8 surrounding voxel centers:

```
f = (p - bmin) / voxel - 0.5           // continuous index
i₀ = floor(fₓ),  u = fₓ - i₀          // integer part + fractional
φ = lerp_3d(φ₀₀₀, φ₁₀₀, ..., φ₁₁₁, u, v, w)
```

Gradient via finite differences of the interpolated values:

```
∂φ/∂x = ((φ₁₀₀ - φ₀₀₀)(1-v)(1-w) + (φ₁₁₀ - φ₀₁₀)v(1-w)
       + (φ₁₀₁ - φ₀₀₁)(1-v)w + (φ₁₁₁ - φ₀₁₁)vw) / hₓ
```

The raw gradient norm `q = |∇φ|` is critical for correct contact force computation — the normal pressure measure is `p(g) · q · dA`, not simply `p(g) · dA`.

### Tricubic

Hermite cubic interpolation using `φ₀` and `n₀` (as gradient) at 8 cell corners. Provides C¹ continuity across cell boundaries, smoother than trilinear.

### Contact-Aware

Precomputes multiple surface-feature candidates per voxel during build. At query time, evaluates each candidate's linear extrapolation and selects the smallest-absolute-gap candidate. Enables multi-contact resolution in regions where multiple surface features are near (e.g., thin features, edges).

## Contact Integration (Python Prototype)

The Python prototype (see `docs/参考意见v1.md`) validates the contact integration pipeline:

1. **Query**: For each integration point on body A's surface, query body B's SDF → `g`, `∇φ`, `|∇φ|`, `n`
2. **Normal pressure**: `p(g) = k_n · max(ε - g, 0)`; effective pressure `p̂ = p · |∇φ|`
3. **Normal traction**: `f_n = p(g) · ∇φ` (uses **raw gradient**, not normalized normal)
4. **Sliding friction**: `f_t = -μ · p̂ · u_t / |u_t|`
5. **Integration**: `F = Σ w_q · (f_n + f_t)`, `τ = Σ w_q · (x - r) × (f_n + f_t)`

Key insight: Friction must be evaluated **per quadrature point** before integration. This preserves the torsional friction moment that would be lost if forces were averaged first.

Last verified against: `plugins/SdfOracle/src/SdfOracle.cpp` (trilinearQuery, firstOrderQuery, etc.), `docs/参考意见v1.md`
