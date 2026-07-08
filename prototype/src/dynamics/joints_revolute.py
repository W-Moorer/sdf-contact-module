import numpy as np


def revolute_joint_residual(joint, model, state):
    frame_i = model.frames[joint.frame_i_id]
    frame_j = model.frames[joint.frame_j_id]

    p_i, A_i = state.world_frame(frame_i, model)
    p_j, A_j = state.world_frame(frame_j, model)

    pos_res = p_i - p_j

    x_i = A_i[:, 0]
    y_i = A_i[:, 1]
    z_j = A_j[:, 2]

    ax_res1 = np.dot(x_i, z_j)
    ax_res2 = np.dot(y_i, z_j)

    return np.concatenate([pos_res, [ax_res1, ax_res2]])


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
