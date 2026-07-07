# Figures

All figures are rendered in Times New Roman font.

## contact_pressure_top.png

Top view of the annular contact patch with effective pressure ($\hat{p}_q = p(g_h)|\nabla\phi_h|$) shown in color. The ring shape (inner radius 0.4, outer radius 0.6) confirms the annular contact patch.

- **Generator**: `prototype/src/visualization.py` — `_plot_top_pressure()`
- **Scene**: `prototype/src/scenes/ring_on_cube_torsion.py` — `run_level0()`
- **Data source**: analytic plane SDF (Level 0), 5120 quadrature points on annular face

## friction_vectors_top.png

Friction traction vectors from top view. Each arrow shows the per-point Coulomb friction traction $f_t = -\mu\hat{p}\,u_t/|u_t|$. The azimuthal (tangential) direction demonstrates that friction forces do not cancel to zero when summed as torques, even though the net force integrates to zero.

- **Generator**: `prototype/src/visualization.py` — `_plot_friction_vectors()`
- **Scene**: `prototype/src/scenes/ring_on_cube_torsion.py` — `run_level0()`
- **Data source**: analytic plane SDF (Level 0), 5120 quadrature points

## torque_density_top.png

Per-point torque density $dT_z = w_q \left[(x_q - r_A) \times f_{t,q}\right]\cdot e_z$. The same-sign (negative) torque contribution across the entire ring proves that torsional friction torque survives integration, even though the tangential force distribution is anti-symmetric.

- **Generator**: `prototype/src/visualization.py` — `_plot_torque_density()`
- **Scene**: `prototype/src/scenes/ring_on_cube_torsion.py` — `run_level0()`
- **Data source**: analytic plane SDF (Level 0), 5120 quadrature points

## torque_error_vs_mesh_resolution.png

Log-log convergence plot of relative torque error vs number of quadrature points. Five OBJ mesh variants (coarse to very_dense) are evaluated against the analytic plane SDF. Error decreases monotonically with mesh refinement.

- **Generator**: `prototype/src/visualization.py` — `save_convergence_plot()`
- **Scene**: `prototype/src/scenes/ring_on_cube_torsion.py` — `run_level2_from_obj()`
- **Input meshes**: `models/hollow_cylinder_{coarse,medium,fine,dense,very_dense}.obj`
- **Generating script**: `scripts/gen_meshes.py`
- **Data source**: analytic plane SDF (Level 0), 5 mesh densities

## Regenerating

Run the Level 0/1 validation suite to regenerate all figures:

```bash
cd prototype
python run_level0.py
```

Output figures are saved under `prototype/outputs/`. Copy the desired ones to `assets/figures/` after regeneration.
