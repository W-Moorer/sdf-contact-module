import numpy as np

class RigidBody:
    def __init__(self, r=np.zeros(3), R=np.eye(3),
                 v=np.zeros(3), omega=np.zeros(3)):
        self.r = np.asarray(r, dtype=float)
        self.R = np.asarray(R, dtype=float)
        self.v = np.asarray(v, dtype=float)
        self.omega = np.asarray(omega, dtype=float)

    def world_from_local(self, X):
        return self.R @ X + self.r

    def local_from_world(self, x):
        return self.R.T @ (x - self.r)

    def velocity_at(self, x):
        return self.v + np.cross(self.omega, x - self.r)
