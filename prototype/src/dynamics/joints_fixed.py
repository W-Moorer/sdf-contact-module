import numpy as np


def fixed_joint_residual(joint, model, state):
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]

    p_i, A_i = state.world_frame(frame_i, model)
    p_j, A_j = state.world_frame(frame_j, model)

    pos_res = p_i - p_j
    R_err = A_j.T @ A_i
    orientation_res = _log_map(R_err)

    return np.concatenate([pos_res, orientation_res])


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
