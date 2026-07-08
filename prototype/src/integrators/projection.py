import numpy as np
from ..dynamics.constraints import _all_residuals, num_joint_constraints, num_drive_constraints
from ..dynamics.mass_matrix import build_mass_matrix
from ..dynamics.utils import _perturb_state
from ..dynamics.joints_fixed import fixed_joint_jacobian
from ..dynamics.joints_revolute import revolute_joint_jacobian
from ..solvers.kkt import solve_mass


def _build_jacobian_full(model, state, use_analytic=True):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return np.zeros((nc, 6 * nb_m))

    if use_analytic:
        J = np.zeros((nc, 6 * nb_m))
        row = 0
        for joint in model.joints.values():
            if joint.type.name == 'FIXED':
                Jj = fixed_joint_jacobian(joint, model, state)
            elif joint.type.name == 'REVOLUTE':
                Jj = revolute_joint_jacobian(joint, model, state)
            else:
                continue
            nr = Jj.shape[0]
            J[row:row+nr, :] = Jj
            row += nr

        nd = num_drive_constraints(model)
        if nd > 0:
            from ..dynamics.drives import motion_drive_jacobian
            Jd = np.zeros((nd, 6 * nb_m))
            r = 0
            for drive in model.drives.values():
                Jd[r:r+1] = motion_drive_jacobian(drive, model, state)
                r += 1
            J[row:row+nd, :] = Jd

        return J
    else:
        eps = 1e-8
        J = np.zeros((nc, 6 * nb_m))
        col = 0
        for body in model.movable_bodies:
            idx = state.body_idx(body.id)
            r_save = state.r[idx].copy()
            R_save = state.R[idx].copy()
            for dof in range(6):
                _perturb_state(state, idx, dof, eps)
                r_p = _all_residuals(model, state)
                state.r[idx] = r_save.copy()
                state.R[idx] = R_save.copy()
                _perturb_state(state, idx, dof, -eps)
                r_m = _all_residuals(model, state)
                state.r[idx] = r_save.copy()
                state.R[idx] = R_save.copy()
                J[:, col] = (r_p - r_m) / (2 * eps)
                col += 1
        return J


def _pseudo_inverse_solve(W, rhs, tol=1e-9):
    U, S, Vt = np.linalg.svd(W, full_matrices=False)
    smax = S[0] if len(S) > 0 else 0.0
    keep = S > tol * max(smax, 1.0)
    return Vt.T[:, keep] @ ((U[:, keep].T @ rhs) / S[keep])


def position_projection(model, state, tol=1e-10, max_iter=20,
                        J_precomputed=None, M_precomputed=None, M_inv_precomputed=None):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return

    for iteration in range(max_iter):
        Phi = _all_residuals(model, state)
        err = np.max(np.abs(Phi))
        if err < tol:
            break

        if iteration == 0 and J_precomputed is not None:
            J = J_precomputed
        else:
            J = _build_jacobian_full(model, state)

        if iteration == 0 and M_precomputed is not None:
            M = M_precomputed
            M_inv = M_inv_precomputed if M_inv_precomputed is not None else np.linalg.inv(M)
        else:
            M = build_mass_matrix(model, state).toarray()
            M_inv = np.linalg.inv(M)

        W = J @ M_inv @ J.T
        lam_p = _pseudo_inverse_solve(W, -Phi)
        delta_V = M_inv @ J.T @ lam_p

        tmp_idx = 0
        for body in model.movable_bodies:
            idx = state.body_idx(body.id)
            state.r[idx] += delta_V[6*tmp_idx:6*tmp_idx+3]
            dtheta = delta_V[6*tmp_idx+3:6*tmp_idx+6]
            dtheta_norm = np.linalg.norm(dtheta)
            if dtheta_norm > 1e-30:
                axis = dtheta / dtheta_norm
                c, s = np.cos(dtheta_norm), np.sin(dtheta_norm)
                dR = c*np.eye(3) + s*np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]) + (1-c)*np.outer(axis,axis)
                state.R[idx] = dR @ state.R[idx]
            tmp_idx += 1


def velocity_projection(model, state,
                        J_precomputed=None, M_precomputed=None, M_inv_precomputed=None):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return

    J = J_precomputed if J_precomputed is not None else _build_jacobian_full(model, state)
    V = state.pack_V()
    M = M_precomputed if M_precomputed is not None else build_mass_matrix(model, state).toarray()
    M_inv = M_inv_precomputed if M_inv_precomputed is not None else np.linalg.inv(M)
    lam_v = _pseudo_inverse_solve(J @ M_inv @ J.T, -(J @ V))
    state.unpack_V(V + M_inv @ J.T @ lam_v)
