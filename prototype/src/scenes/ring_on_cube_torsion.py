import numpy as np
import json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.sdf_grid import AnalyticPlaneSDF, TrilinearSDFGrid
from src.mesh_quadrature import annular_quadrature, load_obj_triangles, triangle_centroid_quadrature, select_bottom_face
from src.rigid_body import RigidBody
from src.contact_integrator import integrate_contact
from src.analytic_benchmarks import analytic_torsion_torque, analytic_equivalent_radius
from src.visualization import save_all_figures, save_convergence_plot, save_sdf_convergence_plot

# ---- shared parameters ----
R_i, R_o = 0.4, 0.6
penetration = 0.002
mu, k_n, omega = 0.5, 1e6, 1.0
params = {'k_n': k_n, 'mu': mu}
models_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'models')

def _body_A():
    return RigidBody(r=[0, 0, -penetration], omega=[0, 0, omega])

def _body_B():
    return RigidBody()

def _run(X_q, w_q, sdf_B, label, extra=None):
    result = integrate_contact(X_q, w_q, _body_A(), _body_B(), sdf_B, params)
    p_uniform = k_n * penetration
    N_analytic = p_uniform * np.pi * (R_o**2 - R_i**2)
    T_z_analytic, R_spin_analytic = analytic_torsion_torque(N_analytic, mu, R_o, R_i, omega)
    R_spin_numeric = (np.sum(w_q * result['p_hat'] * np.linalg.norm(X_q[:, :2], axis=1))
                      / np.sum(w_q * result['p_hat']))
    diag = {
        'label': label,
        'N_numeric': float(result['N_numeric']),
        'N_analytic': float(N_analytic),
        'Tz_numeric': float(result['T_z_numeric']),
        'Tz_analytic': float(T_z_analytic),
        'relative_error_Tz': float(abs(result['T_z_numeric'] - T_z_analytic) / abs(T_z_analytic) if abs(T_z_analytic) > 0 else 0.0),
        'Rspin_numeric': float(R_spin_numeric),
        'Rspin_analytic': float(R_spin_analytic),
        'Ft_residual': float(result['F_t_residual']),
        'Ft_residual_normalized': float(result['F_t_residual'] / (mu * result['N_numeric']) if result['N_numeric'] > 0 else 0.0),
        'force_balance_error': float(np.linalg.norm(result['F_A'] + result['F_B'])),
        'num_quadrature_points': int(len(X_q)),
        'num_active_points': int(np.sum(result['p_hat'] > 0)),
        'R_i': R_i, 'R_o': R_o, 'mu': mu,
    }
    if extra:
        diag.update(extra)
    return diag, result

# =========================================================================
# Level 0 — analytic plane SDF
# =========================================================================
def run_level0(output_dir):
    X_q, w_q, _ = annular_quadrature(R_i, R_o, n_radial=40, n_azimuthal=128)
    sdf_B = AnalyticPlaneSDF(z_top=0.0)
    diag, result = _run(X_q, w_q, sdf_B, "analytic_plane")
    _save(output_dir, 'level0_analytic', diag, result, X_q, w_q)
    return diag

# =========================================================================
# Level 2 — OBJ mesh with analytic plane SDF (mesh convergence)
# =========================================================================
def run_level2_obj(output_dir, obj_path, label):
    verts, tris = load_obj_triangles(obj_path)
    X_q_all, w_q_all, tri_normals, _ = triangle_centroid_quadrature(verts, tris)
    mask = select_bottom_face(X_q_all, tri_normals)
    if mask.sum() == 0:
        print(f"  WARNING: no bottom faces in {label}")
        return None
    X_q, w_q = X_q_all[mask], w_q_all[mask]
    sdf_B = AnalyticPlaneSDF(z_top=0.0)
    diag, result = _run(X_q, w_q, sdf_B, label, {'obj_faces': len(tris), 'bottom_faces': int(mask.sum())})
    _save(output_dir, label, diag, result, X_q, w_q)
    return diag

# =========================================================================
# Level 1 — analytic quadrature with trilinear SDF grid
# =========================================================================
def run_level1_sdf(output_dir, sdf_path, res_label):
    X_q, w_q, _ = annular_quadrature(R_i, R_o, n_radial=40, n_azimuthal=128)
    sdf_B = TrilinearSDFGrid(sdf_path)
    diag, result = _run(X_q, w_q, sdf_B, f"sdf_{res_label}")
    _save(output_dir, f"level1_{res_label}", diag, result, X_q, w_q)
    return diag

# =========================================================================
# Level 1b — OBJ mesh with trilinear SDF grid (combined convergence)
# =========================================================================
def run_level1b_obj(output_dir, obj_path, sdf_path, label, res_label):
    verts, tris = load_obj_triangles(obj_path)
    X_q_all, w_q_all, tri_normals, _ = triangle_centroid_quadrature(verts, tris)
    mask = select_bottom_face(X_q_all, tri_normals)
    if mask.sum() == 0:
        return None
    X_q, w_q = X_q_all[mask], w_q_all[mask]
    sdf_B = TrilinearSDFGrid(sdf_path)
    l = f"{label}_sdf{res_label}"
    diag, result = _run(X_q, w_q, sdf_B, l, {'obj_faces': len(tris), 'bottom_faces': int(mask.sum()), 'sdf_res': res_label})
    _save(output_dir, l, diag, result, X_q, w_q)
    return diag

def _save(base_dir, label, diag, result, X_q, w_q):
    out = os.path.join(base_dir, label)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'diagnostics.json'), 'w') as f:
        json.dump(diag, f, indent=2)
    save_params = {**params, 'R_i': R_i, 'R_o': R_o, 'mu': mu}
    save_all_figures(X_q, w_q, result, save_params, os.path.join(out, 'figures'))
    print(f"  {label:40s}  Tz_err={diag['relative_error_Tz']:.2e}  "
          f"Tz={diag['Tz_numeric']:.4f}  Rspin={diag['Rspin_numeric']:.6f}  "
          f"pts={diag['num_quadrature_points']}")

def run_all(base_dir='outputs'):
    print("=" * 60)
    print("Level 0: analytic plane SDF (reference)")
    print("=" * 60)
    run_level0(os.path.join(base_dir, 'ring_cube_level1'))

    print("\n" + "=" * 60)
    print("Level 2: OBJ mesh + analytic plane SDF (mesh convergence)")
    print("=" * 60)
    obj_variants = ['coarse', 'medium', 'fine', 'dense', 'very_dense']
    mesh_results = []
    for v in obj_variants:
        p = os.path.join(models_dir, f"hollow_cylinder_{v}.obj")
        if os.path.exists(p):
            d = run_level2_obj(base_dir, p, f"mesh_{v}")
            if d: mesh_results.append(d)

    print("\n" + "=" * 60)
    print("Level 1: analytic quadrature + trilinear SDF (SDF convergence)")
    print("=" * 60)
    sdf_results = []
    for res in [64, 128, 256]:
        p = os.path.join(models_dir, f"cube_res{res}.sdf")
        if os.path.exists(p):
            d = run_level1_sdf(base_dir, p, f"res{res}")
            if d: sdf_results.append(d)

    print("\n" + "=" * 60)
    print("Level 1b: OBJ mesh + trilinear SDF (combined convergence)")
    print("=" * 60)
    combined = []
    for v in ['coarse', 'medium', 'fine', 'dense']:
        obj_p = os.path.join(models_dir, f"hollow_cylinder_{v}.obj")
        if not os.path.exists(obj_p): continue
        for res in [64, 128, 256]:
            sdf_p = os.path.join(models_dir, f"cube_res{res}.sdf")
            if os.path.exists(sdf_p):
                d = run_level1b_obj(base_dir, obj_p, sdf_p, v, res)
                if d: combined.append(d)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'Label':>40s}  {'Tz_err':>10s}  {'Tz':>10s}  {'Rspin':>10s}")
    for r in mesh_results + sdf_results:
        print(f"{r['label']:>40s}  {r['relative_error_Tz']:>10.2e}  {r['Tz_numeric']:>10.4f}  {r['Rspin_numeric']:>10.6f}")

    print("\nGenerating convergence plots...")
    base = os.path.join(base_dir, 'ring_cube_level1')
    os.makedirs(os.path.join(base, 'figures'), exist_ok=True)

    if sdf_results:
        save_sdf_convergence_plot(sdf_results, base)
    if mesh_results:
        save_convergence_plot(mesh_results, base, 'mesh')
    if combined:
        save_convergence_plot(combined, base, 'combined')

    return mesh_results, sdf_results, combined

if __name__ == '__main__':
    run_all()
