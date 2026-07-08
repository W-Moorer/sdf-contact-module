import numpy as np
from .joints_fixed import fixed_joint_residual
from .joints_revolute import revolute_joint_residual, revolute_joint_coordinate
from .drives import motion_drive_residual
from .utils import _perturb_state


def eval_constraints(model, state):
    residuals = []
    for joint in model.joints.values():
        if joint.type.name == 'FIXED':
            r = fixed_joint_residual(joint, model, state)
        elif joint.type.name == 'REVOLUTE':
            r = revolute_joint_residual(joint, model, state)
        else:
            raise NotImplementedError(f"Joint type {joint.type} not implemented")
        residuals.append(r)
    if len(residuals) == 0:
        return np.zeros(0)
    return np.concatenate(residuals)


def num_joint_constraints(model):
    n = 0
    for joint in model.joints.values():
        if joint.type.name == 'FIXED':
            n += 6
        elif joint.type.name == 'REVOLUTE':
            n += 5
    return n


def eval_constraint_jacobian(model, state):
    nb_m = model.num_movable
    nc = num_joint_constraints(model)
    if nc == 0 or nb_m == 0:
        return np.zeros((nc, 6 * nb_m))

    J = np.zeros((nc, 6 * nb_m))
    row = 0
    for joint in model.joints.values():
        if joint.type.name == 'FIXED':
            from .joints_fixed import fixed_joint_jacobian
            Jj = fixed_joint_jacobian(joint, model, state)
            nr = Jj.shape[0]
            J[row:row+nr, :] = Jj
            row += nr
        elif joint.type.name == 'REVOLUTE':
            from .joints_revolute import revolute_joint_jacobian
            Jj = revolute_joint_jacobian(joint, model, state)
            nr = Jj.shape[0]
            J[row:row+nr, :] = Jj
            row += nr
        else:
            raise NotImplementedError(f"Analytic Jacobian for {joint.type}")

    return J


def eval_drive_constraints(model, state):
    residuals = []
    for drive in model.drives.values():
        r = motion_drive_residual(drive, model, state)
        if r is not None:
            residuals.append(r)
    if len(residuals) == 0:
        return np.zeros(0)
    return np.concatenate(residuals)


def num_drive_constraints(model):
    return len(model.drives)


def compute_bc(model, state, dt, alpha=0.1, beta=0.1):
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0:
        return np.zeros(nc)

    eps = 1e-8
    J_tilde = np.zeros((nc, 6 * model.num_movable))

    col = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        r_save = state.r[idx].copy()
        R_save = state.R[idx].copy()

        for dof in range(6):
            _perturb_state(state, idx, dof, eps)
            Phi_plus = _all_residuals(model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()
            _perturb_state(state, idx, dof, -eps)
            Phi_minus = _all_residuals(model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()

            J_tilde[:, col] = (Phi_plus - Phi_minus) / (2 * eps)
            col += 1

    Phi = _all_residuals(model, state)

    V = state.pack_V()
    Phi_dot = J_tilde @ V

    nd_joint = num_joint_constraints(model)
    for i, drive in enumerate(model.drives.values()):
        func = model.functions.get(drive.function_id)
        if func is not None:
            Phi_dot[nd_joint + i] += -func.first_derivative(state.t)

    if dt > 0:
        b_stab = 2 * alpha * Phi_dot + beta**2 * Phi
    else:
        b_stab = np.zeros(nc)

    b_c = -Phi_dot - b_stab
    return b_c


def compute_bc_from_jacobian(J, model, state, dt, alpha=0.1, beta=0.1):
    nc = J.shape[0]
    if nc == 0:
        return np.zeros(0)

    Phi = _all_residuals(model, state)
    V = state.pack_V()
    Phi_dot = J @ V

    nd_joint = num_joint_constraints(model)
    for i, drive in enumerate(model.drives.values()):
        func = model.functions.get(drive.function_id)
        if func is not None:
            Phi_dot[nd_joint + i] += -func.first_derivative(state.t)

    if dt > 0:
        b_stab = 2 * alpha * Phi_dot + beta**2 * Phi
    else:
        b_stab = np.zeros(nc)

    b_c = -Phi_dot - b_stab
    return b_c


def _all_residuals(model, state):
    joint_res = eval_constraints(model, state)
    drive_res = eval_drive_constraints(model, state)
    if len(joint_res) == 0 and len(drive_res) == 0:
        return np.zeros(0)
    return np.concatenate([joint_res, drive_res])
