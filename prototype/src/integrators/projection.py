import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from ..dynamics.constraints import _all_residuals, num_joint_constraints, num_drive_constraints
from ..dynamics.mass_matrix import build_mass_matrix
from ..dynamics.utils import _perturb_state
from ..solvers.kkt import solve_mass


def _build_jacobian_full(model, state):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return np.zeros((nc, 6 * nb_m))

    eps = 1e-8
    J = np.zeros((nc, 6 * nb_m))

    col = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        r_save = state.r[idx].copy()
        R_save = state.R[idx].copy()

        for dof in range(6):
            _perturb_state(state, idx, dof, eps)
            r_plus = _all_residuals(model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()
            _perturb_state(state, idx, dof, -eps)
            r_minus = _all_residuals(model, state)

            state.r[idx] = r_save.copy()
            state.R[idx] = R_save.copy()

            J[:, col] = (r_plus - r_minus) / (2 * eps)
            col += 1

    return J


def _pseudo_inverse_solve(W, rhs, tol=1e-9):
    """Solve W @ x = rhs using SVD pseudo-inverse, handling rank deficiency."""
    U, S, Vt = np.linalg.svd(W, full_matrices=False)
    smax = S[0] if len(S) > 0 else 0.0
    keep = S > tol * max(smax, 1.0)
    return Vt.T[:, keep] @ ((U[:, keep].T @ rhs) / S[keep])


def position_projection(model, state, tol=1e-10, max_iter=20):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return

    for iteration in range(max_iter):
        Phi = _all_residuals(model, state)
        err = np.max(np.abs(Phi))
        if err < tol:
            break

        J = _build_jacobian_full(model, state)
        M = build_mass_matrix(model, state).toarray()
        M_inv = np.linalg.inv(M)

        # Use J @ M^{-1} @ J^T with pseudo-inverse
        W = J @ M_inv @ J.T
        rhs = -Phi
        lam_p = _pseudo_inverse_solve(W, rhs)

        # delta_q = M^{-1} @ J^T @ lam_p
        delta_V = M_inv @ J.T @ lam_p

        tmp_idx = 0
        for body in model.movable_bodies:
            idx = state.body_idx(body.id)
            state.r[idx] += delta_V[6*tmp_idx:6*tmp_idx+3]

            dtheta = delta_V[6*tmp_idx+3:6*tmp_idx+6]
            dtheta_norm = np.linalg.norm(dtheta)
            if dtheta_norm > 1e-30:
                axis = dtheta / dtheta_norm
                c = np.cos(dtheta_norm)
                s = np.sin(dtheta_norm)
                dR = (c * np.eye(3)
                      + s * np.array([[0, -axis[2], axis[1]],
                                      [axis[2], 0, -axis[0]],
                                      [-axis[1], axis[0], 0]])
                      + (1 - c) * np.outer(axis, axis))
                state.R[idx] = dR @ state.R[idx]
            tmp_idx += 1


def velocity_projection(model, state, tol=1e-9):
    nb_m = model.num_movable
    nc = num_joint_constraints(model) + num_drive_constraints(model)
    if nc == 0 or nb_m == 0:
        return

    J = _build_jacobian_full(model, state)
    V = state.pack_V()
    rhs = -(J @ V)

    M = build_mass_matrix(model, state).toarray()
    M_inv = np.linalg.inv(M)
    W = J @ M_inv @ J.T

    lam_v = _pseudo_inverse_solve(W, rhs)
    delta_V = M_inv @ J.T @ lam_v
    state.unpack_V(V + delta_V)
