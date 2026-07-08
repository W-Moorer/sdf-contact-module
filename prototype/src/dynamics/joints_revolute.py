import numpy as np


def _skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def _body_col(model, body):
    for i, b in enumerate(model.movable_bodies):
        if b.id == body.id:
            return 6 * i
    return -1


def revolute_joint_residual(joint, model, state):
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]
    p_i, A_i = state.world_frame(frame_i, model)
    p_j, A_j = state.world_frame(frame_j, model)

    x_i = A_i[:, 0]
    y_i = A_i[:, 1]
    z_j = A_j[:, 2]

    return np.concatenate([p_i - p_j, [np.dot(x_i, z_j), np.dot(y_i, z_j)]])


def revolute_joint_jacobian(joint, model, state):
    """Analytic Jacobian of RevoluteJoint: 5 rows x (6*n_bodies) cols.

    Rows:
      0-2: position p_i - p_j
      3:   x_i · z_j = 0
      4:   y_i · z_j = 0
    """
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]
    body_i = model.bodies[frame_i.body_id]
    body_j = model.bodies[frame_j.body_id]

    idx_i = state.body_idx(body_i.id)
    idx_j = state.body_idx(body_j.id)
    nb_m = model.num_movable
    J = np.zeros((5, 6 * nb_m))

    R_i = state.R[idx_i]
    R_j = state.R[idx_j]
    s_i = frame_i.local_position
    s_j = frame_j.local_position
    R_i0 = frame_i.local_orientation
    R_j0 = frame_j.local_orientation

    col_i = _body_col(model, body_i)
    col_j = _body_col(model, body_j)

    # =========== Position constraints (rows 0-2) ===========
    # ∂Φ_p/∂r_i = I, ∂Φ_p/∂θ_i = -R_i·[s_i]×
    # ∂Φ_p/∂r_j = -I, ∂Φ_p/∂θ_j = R_j·[s_j]×
    if col_i >= 0:
        J[0:3, col_i:col_i+3] = np.eye(3)
        J[0:3, col_i+3:col_i+6] = -R_i @ _skew(s_i)
    if col_j >= 0:
        J[0:3, col_j:col_j+3] = -np.eye(3)
        J[0:3, col_j+3:col_j+6] = R_j @ _skew(s_j)

    # =========== Axis alignment constraints (rows 3-4) ===========
    # x_i = A_i[:, 0] = R_i @ R_i0[:, 0],  z_j = A_j[:, 2] = R_j @ R_j0[:, 2]
    # y_i = A_i[:, 1] = R_i @ R_i0[:, 1]
    x_i = R_i @ R_i0[:, 0]
    y_i = R_i @ R_i0[:, 1]
    z_j = R_j @ R_j0[:, 2]

    # ∂(x_i·z_j)/∂θ_i = -(x_i × z_j)   [row 3, cols for body i]
    # ∂(x_i·z_j)/∂θ_j =  (x_i × z_j)   [row 3, cols for body j]
    # ∂(y_i·z_j)/∂θ_i = -(y_i × z_j)   [row 4, cols for body i]
    # ∂(y_i·z_j)/∂θ_j =  (y_i × z_j)   [row 4, cols for body j]
    x_cross_z = np.cross(x_i, z_j)
    y_cross_z = np.cross(y_i, z_j)

    if col_i >= 0:
        J[3, col_i+3:col_i+6] = -x_cross_z
        J[4, col_i+3:col_i+6] = -y_cross_z
    if col_j >= 0:
        J[3, col_j+3:col_j+6] = x_cross_z
        J[4, col_j+3:col_j+6] = y_cross_z

    return J


def revolute_joint_coordinate(joint, model, state):
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]
    _, A_i = state.world_frame(frame_i, model)
    _, A_j = state.world_frame(frame_j, model)

    x_i = A_i[:, 0]
    y_i = A_i[:, 1]
    x_j = A_j[:, 0]
    y_j = A_j[:, 1]

    theta = np.arctan2(np.dot(x_i, y_j), np.dot(x_i, x_j))
    return theta


def revolute_joint_coordinate_jacobian(joint, model, state):
    """Finite-difference Jacobian of the revolute joint coordinate."""
    from .utils import _perturb_state
    eps = 1e-8
    theta0 = revolute_joint_coordinate(joint, model, state)
    nb_m = model.num_movable
    J = np.zeros(6 * nb_m)

    col = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        r_save = state.r[idx].copy()
        R_save = state.R[idx].copy()
        for dof in range(6):
            _perturb_state(state, idx, dof, eps)
            theta_p = revolute_joint_coordinate(joint, model, state)
            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()
            _perturb_state(state, idx, dof, -eps)
            theta_m = revolute_joint_coordinate(joint, model, state)
            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()
            J[col] = (theta_p - theta_m) / (2 * eps)
            col += 1

    return J
