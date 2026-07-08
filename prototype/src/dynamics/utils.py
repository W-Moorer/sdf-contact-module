import numpy as np


def _exp_map(axis_angle):
    theta = np.linalg.norm(axis_angle)
    if theta < 1e-30:
        return np.eye(3)
    axis = axis_angle / theta
    c = np.cos(theta)
    s = np.sin(theta)
    return (c * np.eye(3)
            + s * np.array([[0, -axis[2], axis[1]],
                            [axis[2], 0, -axis[0]],
                            [-axis[1], axis[0], 0]])
            + (1 - c) * np.outer(axis, axis))


def _perturb_state(state, body_idx, dof, eps):
    if dof < 3:
        state.r[body_idx, dof] += eps
    else:
        delta = np.zeros(3)
        delta[dof - 3] = eps
        dR = _exp_map(delta)
        state.R[body_idx] = dR @ state.R[body_idx]
