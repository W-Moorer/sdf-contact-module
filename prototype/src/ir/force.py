import numpy as np


class ForceElement:
    def __init__(self, name, force_id, frame_i_id, frame_j_id,
                 magnitude_function=None, force_type='axial'):
        self.name = name
        self.id = force_id
        self.frame_i_id = frame_i_id
        self.frame_j_id = frame_j_id
        self.magnitude_function = magnitude_function
        self.force_type = force_type

    def evaluate(self, state, model):
        frame_i = model.frames[self.frame_i_id]
        frame_j = model.frames[self.frame_j_id]

        p_i, _ = state.world_frame(frame_i, model)
        p_j, _ = state.world_frame(frame_j, model)

        body_i = model.bodies[frame_i.body_id]
        body_j = model.bodies[frame_j.body_id]
        idx_i = state.body_idx(body_i.id)
        idx_j = state.body_idx(body_j.id)

        f_mag = self.magnitude_function.value(state.t) if self.magnitude_function else 0.0

        if self.force_type == 'axial':
            d = p_j - p_i
            dist = np.linalg.norm(d)
            if dist > 1e-30:
                direction = d / dist
            else:
                direction = np.zeros(3)
            F_i = f_mag * direction
            F_j = -F_i
        else:
            F_i = np.zeros(3)
            F_j = np.zeros(3)

        tau_i = np.cross(p_i - state.r[idx_i], F_i)
        tau_j = np.cross(p_j - state.r[idx_j], F_j)

        return F_i, tau_i, F_j, tau_j


class TorqueElement:
    def __init__(self, name, torque_id, body_id, axis_local=None,
                 magnitude_function=None):
        self.name = name
        self.id = torque_id
        self.body_id = body_id
        self.axis_local = np.array([0, 0, 1]) if axis_local is None else np.array(axis_local, dtype=float)
        self.magnitude_function = magnitude_function

    def evaluate(self, state, model):
        body = model.bodies[self.body_id]
        idx = state.body_idx(body.id)
        t_mag = self.magnitude_function.value(state.t) if self.magnitude_function else 0.0
        tau = t_mag * (state.R[idx] @ self.axis_local)
        return tau
