import numpy as np
import json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.sdf_grid import AnalyticPlaneSDF
from src.mesh_quadrature import annular_quadrature
from src.rigid_body import RigidBody
from src.contact_integrator import integrate_contact
from src.analytic_benchmarks import analytic_torsion_torque, analytic_equivalent_radius
from src.visualization import save_all_figures

def run_ring_on_cube_torsion(output_dir='outputs/ring_cube_level0'):
    R_i = 0.4
    R_o = 0.6
    penetration = 0.002
    mu = 0.5
    k_n = 1e6
    omega = 1.0
    epsilon = 0.0
    regularizer = 0.0

    sdf_B = AnalyticPlaneSDF(z_top=0.0)

    X_q, w_q, tri_id = annular_quadrature(R_i, R_o, n_radial=20, n_azimuthal=64)

    body_A = RigidBody(
        r=np.array([0.0, 0.0, -penetration]),
        R=np.eye(3),
        v=np.zeros(3),
        omega=np.array([0.0, 0.0, omega])
    )
    body_B = RigidBody(
        r=np.zeros(3),
        R=np.eye(3),
        v=np.zeros(3),
        omega=np.zeros(3)
    )

    params = {'k_n': k_n, 'mu': mu, 'epsilon': epsilon, 'regularizer': regularizer}
    result = integrate_contact(X_q, w_q, body_A, body_B, sdf_B, params)

    N_numeric = result['N_numeric']
    T_z_numeric = result['T_z_numeric']
    F_t_residual = result['F_t_residual']

    p_uniform = k_n * penetration
    N_analytic = p_uniform * np.pi * (R_o**2 - R_i**2)
    T_z_analytic, R_spin_analytic = analytic_torsion_torque(N_analytic, mu, R_o, R_i, omega)
    R_spin_numeric = (np.sum(w_q * result['p_hat'] * np.linalg.norm(X_q[:, :2], axis=1))
                      / np.sum(w_q * result['p_hat']))

    relative_error = abs(T_z_numeric - T_z_analytic) / abs(T_z_analytic) if abs(T_z_analytic) > 0 else 0.0

    friction_power = np.sum(w_q * np.sum(result['forces_per_point'] * 0, axis=1))
    force_balance_error = np.linalg.norm(result['F_A'] + result['F_B'])

    diagnostics = {
        'N_numeric': float(N_numeric),
        'N_analytic': float(N_analytic),
        'Tz_numeric': float(T_z_numeric),
        'Tz_analytic': float(T_z_analytic),
        'relative_error_Tz': float(relative_error),
        'Rspin_numeric': float(R_spin_numeric),
        'Rspin_analytic': float(R_spin_analytic),
        'Ft_residual': float(F_t_residual),
        'Ft_residual_normalized': float(F_t_residual / (mu * N_numeric) if N_numeric > 0 else 0.0),
        'friction_power': float(friction_power),
        'force_balance_error': float(force_balance_error),
        'num_quadrature_points': int(len(X_q)),
        'num_active_points': int(np.sum(result['p_hat'] > 0)),
        'R_i': R_i, 'R_o': R_o, 'mu': mu, 'penetration': penetration, 'omega': omega
    }

    print(f"N_numeric          = {N_numeric:.6e}")
    print(f"N_analytic         = {N_analytic:.6e}")
    print(f"Tz_numeric         = {T_z_numeric:.6e}")
    print(f"Tz_analytic        = {T_z_analytic:.6e}")
    print(f"relative_error_Tz  = {relative_error:.6e}")
    print(f"Rspin_numeric      = {R_spin_numeric:.6f}")
    print(f"Rspin_analytic     = {R_spin_analytic:.6f}")
    print(f"Ft_residual_norm   = {diagnostics['Ft_residual_normalized']:.6e}")
    print(f"force_balance_err  = {force_balance_error:.6e}")

    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'diagnostics.json'), 'w') as f:
        json.dump(diagnostics, f, indent=2)

    save_params = {**params, 'R_i': R_i, 'R_o': R_o, 'mu': mu}
    save_all_figures(X_q, w_q, result, save_params,
                     os.path.join(output_dir, 'figures'))

    return diagnostics

if __name__ == '__main__':
    run_ring_on_cube_torsion()
