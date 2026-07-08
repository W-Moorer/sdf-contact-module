import numpy as np


class RigidBody:
    def __init__(self, name, body_id, mass=1.0, inertia_body=None,
                 com_local=None, initial_r=None, initial_R=None,
                 initial_v=None, initial_omega=None, fixed=False):
        self.name = name
        self.id = body_id
        self.mass = mass
        self.inertia_body = np.eye(3) if inertia_body is None else np.asarray(inertia_body, dtype=float)
        self.com_local = np.zeros(3) if com_local is None else np.asarray(com_local, dtype=float)
        self.initial_r = np.zeros(3) if initial_r is None else np.asarray(initial_r, dtype=float)
        self.initial_R = np.eye(3) if initial_R is None else np.asarray(initial_R, dtype=float)
        self.initial_v = np.zeros(3) if initial_v is None else np.asarray(initial_v, dtype=float)
        self.initial_omega = np.zeros(3) if initial_omega is None else np.asarray(initial_omega, dtype=float)
        self.fixed = fixed
