import numpy as np


class State:
    def __init__(self, model):
        nb = model.num_bodies
        body_list = list(model.bodies.values())
        self.r = np.zeros((nb, 3))
        self.R = np.zeros((nb, 3, 3))
        self.v = np.zeros((nb, 3))
        self.omega = np.zeros((nb, 3))
        self.t = 0.0

        self._body_index = {}
        for i, b in enumerate(body_list):
            self._body_index[b.id] = i
            self.r[i] = b.initial_r
            self.R[i] = b.initial_R
            self.v[i] = b.initial_v
            self.omega[i] = b.initial_omega

        self._movable_ids = [b.id for b in model.movable_bodies]
        self._nb_all = nb

    def body_idx(self, body_id):
        return self._body_index[body_id]

    def pack_V(self):
        nb_m = len(self._movable_ids)
        V = np.zeros(6 * nb_m)
        for i, bid in enumerate(self._movable_ids):
            idx = self._body_index[bid]
            V[6*i:6*i+3] = self.v[idx]
            V[6*i+3:6*i+6] = self.omega[idx]
        return V

    def unpack_V(self, V):
        nb_m = len(self._movable_ids)
        for i, bid in enumerate(self._movable_ids):
            idx = self._body_index[bid]
            self.v[idx] = V[6*i:6*i+3]
            self.omega[idx] = V[6*i+3:6*i+6]

    def copy(self):
        s = State.__new__(State)
        s.r = self.r.copy()
        s.R = self.R.copy()
        s.v = self.v.copy()
        s.omega = self.omega.copy()
        s.t = self.t
        s._body_index = self._body_index
        s._movable_ids = self._movable_ids
        s._nb_all = self._nb_all
        return s

    def world_frame(self, frame, model):
        idx = self.body_idx(frame.body_id)
        p = self.r[idx] + self.R[idx] @ frame.local_position
        A = self.R[idx] @ frame.local_orientation
        return p, A
