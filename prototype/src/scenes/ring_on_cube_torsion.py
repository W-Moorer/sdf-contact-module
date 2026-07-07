import numpy as np
import json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.sdf_grid import AnalyticPlaneSDF
from src.mesh_quadrature import annular_quadrature, load_obj_triangles, triangle_centroid_quadrature, select_bottom_face
from src.rigid_body import RigidBody
from src.contact_integrator import integrate_contact
from src.analytic_benchmarks import analytic_torsion_torque, analytic_equivalent_radius
from src.visualization import save_all_figures

def run_level0_analytic(output_dir='outputs/ring_cube_level0',
                         n_radial=20, n_azimuthal=64):
    R_i, R_o = 0.4, 0.6
    penetration = 0.002
    mu, k_n, omega = 0.5, 1e6, 1.0

    sdf_B = AnalyticPlaneSDF(z_top=0.0)
    X_q, w_q, _ = annular_quadrature(R_i, R_o, n_radial, n_azimuthal)

    body_A = RigidBody(r=[0, 0, -penetration], omega=[0, 0, omega])
    body_B = RigidBody()

    params = {'k_n': k_n, 'mu': mu}
    result = integrate_contact(X_q, w_q, body_A, body_B, sdf_B, params)

    p_uniform = k_n * penetration
    N_analytic = p_uniform * np.pi * (R_o**2 - R_i**2)
    T_z_analytic, R_spin_analytic = analytic_torsion_torque(N_analytic, mu, R_o, R_i, omega)
    R_spin_numeric = (np.sum(w_q * result['p_hat'] * np.linalg.norm(X_q[:, :2], axis=1))
                      / np.sum(w_q * result['p_hat']))

    diag = _make_diagnostics(result, N_analytic, T_z_analytic, R_spin_analytic,
                             R_spin_numeric, mu, f"analytic_quad_{n_radial}x{n_azimuthal}",
                             R_i, R_o, penetration, omega, params)
    _save_outputs(output_dir, diag, X_q, w_q, result, params, R_i, R_o, mu)
    return diag

def run_level2_from_obj(output_dir, obj_path, label):
    R_i, R_o = 0.4, 0.6
    penetration = 0.002
    mu, k_n, omega = 0.5, 1e6, 1.0

    verts, tris = load_obj_triangles(obj_path)
    X_q_all, w_q_all, tri_normals, _ = triangle_centroid_quadrature(verts, tris)

    mask = select_bottom_face(X_q_all, tri_normals)
    if mask.sum() == 0:
        print(f"  WARNING: no bottom faces found in {label}")
        return None

    X_q = X_q_all[mask]
    w_q = w_q_all[mask]

    sdf_B = AnalyticPlaneSDF(z_top=0.0)
    body_A = RigidBody(r=[0, 0, -penetration], omega=[0, 0, omega])
    body_B = RigidBody()

    params = {'k_n': k_n, 'mu': mu}
    result = integrate_contact(X_q, w_q, body_A, body_B, sdf_B, params)

    p_uniform = k_n * penetration
    N_analytic = p_uniform * np.pi * (R_o**2 - R_i**2)
    T_z_analytic, R_spin_analytic = analytic_torsion_torque(N_analytic, mu, R_o, R_i, omega)
    R_spin_numeric = (np.sum(w_q * result['p_hat'] * np.linalg.norm(X_q[:, :2], axis=1))
                      / np.sum(w_q * result['p_hat']))

    diag = _make_diagnostics(result, N_analytic, T_z_analytic, R_spin_analytic,
                             R_spin_numeric, mu, label,
                             R_i, R_o, penetration, omega, params)
    diag['obj_faces'] = len(tris)
    diag['bottom_faces'] = int(mask.sum())

    _save_outputs(os.path.join(output_dir, label), diag,
                  X_q, w_q, result, params, R_i, R_o, mu)
    return diag

def _make_diagnostics(result, N_analytic, T_z_analytic, R_spin_analytic,
                      R_spin_numeric, mu, label, R_i, R_o, penetration, omega, params):
    T_z_numeric = result['T_z_numeric']
    N_numeric = result['N_numeric']
    relative_error = abs(T_z_numeric - T_z_analytic) / abs(T_z_analytic) if abs(T_z_analytic) > 0 else 0.0
    force_balance_error = np.linalg.norm(result['F_A'] + result['F_B'])

    return {
        'label': label,
        'N_numeric': float(N_numeric),
        'N_analytic': float(N_analytic),
        'Tz_numeric': float(T_z_numeric),
        'Tz_analytic': float(T_z_analytic),
        'relative_error_Tz': float(relative_error),
        'Rspin_numeric': float(R_spin_numeric),
        'Rspin_analytic': float(R_spin_analytic),
        'Ft_residual': float(result['F_t_residual']),
        'Ft_residual_normalized': float(result['F_t_residual'] / (mu * N_numeric) if N_numeric > 0 else 0.0),
        'force_balance_error': float(force_balance_error),
        'num_quadrature_points': int(len(result['gaps'])),
        'num_active_points': int(np.sum(result['p_hat'] > 0)),
        'R_i': R_i, 'R_o': R_o, 'mu': mu, 'penetration': penetration, 'omega': omega,
    }

def _save_outputs(output_dir, diag, X_q, w_q, result, params, R_i, R_o, mu):
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'diagnostics.json'), 'w') as f:
        json.dump(diag, f, indent=2)
    save_params = {**params, 'R_i': R_i, 'R_o': R_o, 'mu': mu}
    save_all_figures(X_q, w_q, result, save_params,
                     os.path.join(output_dir, 'figures'))
    print(f"  {diag['label']:30s}  Tz_err={diag['relative_error_Tz']:.2e}  "
          f"N={diag['N_numeric']:.2f}  Tz={diag['Tz_numeric']:.4f}  "
          f"pts={diag['num_quadrature_points']}")

def run_all_comparisons(base_dir='outputs'):
    print("=" * 60)
    print("Level 0: analytic plane SDF")
    print("=" * 60)

    print("\n--- Analytic quadrature (high-res reference) ---")
    run_level0_analytic(os.path.join(base_dir, 'ring_cube_level0'),
                        n_radial=40, n_azimuthal=128)

    print("\n--- OBJ mesh variants (bottom face only) ---")
    models_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'models')
    variants = ['coarse', 'medium', 'fine', 'dense', 'very_dense']
    results = []
    for v in variants:
        obj_path = os.path.join(models_dir, f"hollow_cylinder_{v}.obj")
        if os.path.exists(obj_path):
            diag = run_level2_from_obj(base_dir, obj_path, v)
            if diag:
                results.append(diag)

    print("\n--- Summary ---")
    print(f"{'Variant':>15s}  {'Faces':>6s}  {'Pts':>6s}  {'Tz_err':>10s}  {'Tz_num':>10s}  {'Rspin_err':>10s}")
    for r in results:
        rsp_err = abs(r['Rspin_numeric'] - r['Rspin_analytic']) / r['Rspin_analytic']
        print(f"{r['label']:>15s}  {r.get('obj_faces',0):>6d}  {r['num_quadrature_points']:>6d}  "
              f"{r['relative_error_Tz']:>10.2e}  {r['Tz_numeric']:>10.4f}  {rsp_err:>10.2e}")

    return results

if __name__ == '__main__':
    run_all_comparisons()
