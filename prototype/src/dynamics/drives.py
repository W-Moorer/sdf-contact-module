import numpy as np
from .joints_revolute import revolute_joint_coordinate
from .utils import _perturb_state


def motion_drive_residual(drive, model, state):
    joint = model.joints.get(drive.target_joint_id)
    if joint is None:
        return None

    if joint.type.name == 'REVOLUTE':
        theta = revolute_joint_coordinate(joint, model, state)
    else:
        return None

    func = model.functions.get(drive.function_id)
    if func is None:
        return np.array([theta])

    theta_d = func.value(state.t)
    return np.array([theta - theta_d])


def motion_drive_jacobian(drive, model, state):
    eps = 1e-8
    nb_m = model.num_movable
    J_drive = np.zeros((1, 6 * nb_m))

    col = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        r_save = state.r[idx].copy()
        R_save = state.R[idx].copy()

        for dof in range(6):
            _perturb_state(state, idx, dof, eps)
            r_plus = motion_drive_residual(drive, model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()
            _perturb_state(state, idx, dof, -eps)
            r_minus = motion_drive_residual(drive, model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()

            J_drive[0, col] = (r_plus[0] - r_minus[0]) / (2 * eps)
            col += 1

    return J_drive
