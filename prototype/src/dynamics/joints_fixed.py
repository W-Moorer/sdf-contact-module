import numpy as np


def _skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def fixed_joint_residual(joint, model, state):
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]
    p_i, A_i = state.world_frame(frame_i, model)
    p_j, A_j = state.world_frame(frame_j, model)

    pos_res = p_i - p_j
    R_err = A_j.T @ A_i
    orientation_res = _log_map(R_err)

    return np.concatenate([pos_res, orientation_res])


def fixed_joint_jacobian(joint, model, state):
    """Analytic Jacobian of FixedJoint: 6 rows x (6*n_bodies) cols."""
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]
    body_i = model.bodies[frame_i.body_id]
    body_j = model.bodies[frame_j.body_id]

    idx_i = state.body_idx(body_i.id)
    idx_j = state.body_idx(body_j.id)
    nb_m = model.num_movable
    J = np.zeros((6, 6 * nb_m))

    R_i = state.R[idx_i]
    R_j = state.R[idx_j]
    s_i = frame_i.local_position
    s_j = frame_j.local_position

    # Position constraints (rows 0-2): Φ_p = p_i - p_j
    # ∂Φ_p/∂r_i = I
    # ∂Φ_p/∂θ_i = -R_i · [s_i]×
    # ∂Φ_p/∂r_j = -I
    # ∂Φ_p/∂θ_j = R_j · [s_j]×
    col_i = _body_col(model, body_i)
    col_j = _body_col(model, body_j)

    if col_i >= 0:
        J[0:3, col_i:col_i+3] = np.eye(3)
        J[0:3, col_i+3:col_i+6] = -R_i @ _skew(s_i)
    if col_j >= 0:
        J[0:3, col_j:col_j+3] = -np.eye(3)
        J[0:3, col_j+3:col_j+6] = R_j @ _skew(s_j)

    # Orientation constraints (rows 3-5): Φ_R = Log(A_j^T @ A_i)
    # A_i = R_i @ R_a0, A_j = R_j @ R_b0
    A_i = R_i @ frame_i.local_orientation
    A_j = R_j @ frame_j.local_orientation
    R_err = A_j.T @ A_i

    # For small orientation errors: dΦ_R/dθ_i ≈ A_j^T, dΦ_R/dθ_j ≈ -A_j^T
    # Using left-trivialized derivative of Log map
    J_ort = _log_map_derivative_inv(R_err)
    AjT = A_j.T

    if col_i >= 0:
        J[3:6, col_i+3:col_i+6] = J_ort @ AjT
    if col_j >= 0:
        J[3:6, col_j+3:col_j+6] = -J_ort @ AjT

    return J


def _body_col(model, body):
    """Return column index in KKT for this body, or -1 if fixed."""
    for i, b in enumerate(model.movable_bodies):
        if b.id == body.id:
            return 6 * i
    return -1


def _log_map(R):
    cos_theta = (np.trace(R) - 1.0) / 2.0
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    if theta < 1e-10:
        return np.zeros(3)
    s = theta / (2.0 * np.sin(theta))
    return s * np.array([R[2, 1] - R[1, 2],
                         R[0, 2] - R[2, 0],
                         R[1, 0] - R[0, 1]])


def _log_map_derivative_inv(R):
    """Left-trivialized derivative inverse of Log map at R.
    
    d(Log(R))/dR = T(R)^{-1} where T(R) = (θ/2·cot(θ/2))·I + 
                    (1 - θ/2·cot(θ/2))·n·n^T - (θ/2)·[n]×
    
    For θ → 0: T^{-1} ≈ I + 0.5·[n]× + ...
    """
    cos_theta = (np.trace(R) - 1.0) / 2.0
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    if theta < 1e-10:
        return np.eye(3)

    n = np.array([R[2, 1] - R[1, 2],
                  R[0, 2] - R[2, 0],
                  R[1, 0] - R[0, 1]])
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-10:
        return np.eye(3)
    n = n / n_norm

    half = theta / 2.0
    cot_half = 1.0 / np.tan(half) if np.tan(half) != 0 else 0.0
    a = half * cot_half
    b = 1.0 - a

    # T = a*I + b*n*n^T - half*[n]×
    # T^{-1} can be computed as:
    # T^{-1} = (1/a)·I - (b/(a*(a+b)))·n·n^T + (1/(2a))·[n]×
    # For θ → 0: a → 1, b → 0, T^{-1} → I + 0.5·[n]×
    # Actually the standard formula for T^{-1}(θ,n) is:
    # T_inv = I + (1/a - 1)·n·n^T + (1/(2a))·[n]×
    T_inv = (1.0 / a) * np.eye(3) \
          + (1.0 / a - 1.0) * np.outer(n, n) \
          + (0.5 / a) * _skew(n)

    return T_inv
