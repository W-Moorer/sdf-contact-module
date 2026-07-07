import numpy as np
from .contact_law import compute_normal_pressure, compute_friction_traction

def integrate_contact(X_q, w_q, body_A, body_B, sdf_B, params):
    k_n = params.get('k_n', 1e6)
    mu = params.get('mu', 0.5)
    epsilon = params.get('epsilon', 0.0)
    regularizer = params.get('regularizer', 0.0)

    N = len(X_q)
    forces_A = np.zeros((N, 3))
    torques_A = np.zeros((N, 3))

    gaps = np.zeros(N)
    p_hats = np.zeros(N)

    for i in range(N):
        x_q = body_A.world_from_local(X_q[i])
        y_q = body_B.local_from_world(x_q)
        g, raw_grad, grad_norm, unit_normal, valid = sdf_B.query(y_q)

        p = compute_normal_pressure(g, k_n, epsilon)
        p_hat = p * grad_norm

        f_n = p * raw_grad

        u_A = body_A.velocity_at(x_q)
        u_B = body_B.velocity_at(x_q)
        u_rel = u_A - u_B
        u_t = u_rel - np.dot(u_rel, unit_normal) * unit_normal

        f_t = compute_friction_traction(u_t, p_hat, mu, regularizer)

        f_q = f_n + f_t
        forces_A[i] = w_q[i] * f_q
        torques_A[i] = w_q[i] * np.cross(x_q - body_A.r, f_q)
        gaps[i] = g
        p_hats[i] = p_hat

    F_A = np.sum(forces_A, axis=0)
    tau_A = np.sum(torques_A, axis=0)
    F_B = -F_A

    total_torque_world = np.sum(torques_A, axis=0)
    tau_B = -total_torque_world + np.cross(body_A.r - body_B.r, F_B)

    result = {
        'F_A': F_A,
        'tau_A': tau_A,
        'F_B': F_B,
        'tau_B': tau_B,
        'gaps': gaps,
        'p_hat': p_hats,
        'forces_per_point': forces_A,
        'torques_per_point': torques_A,
        'T_z_numeric': tau_A[2],
        'N_numeric': np.sum(w_q * p_hats),
        'F_t_residual': np.linalg.norm(F_A[:2]),
    }
    return result
