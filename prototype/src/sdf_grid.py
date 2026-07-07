import numpy as np

class AnalyticPlaneSDF:
    def __init__(self, z_top=0.0):
        self.z_top = z_top

    def query(self, y_local):
        g = y_local[2] - self.z_top
        raw_grad = np.array([0.0, 0.0, 1.0])
        grad_norm = 1.0
        unit_normal = np.array([0.0, 0.0, 1.0])
        valid = True
        return g, raw_grad, grad_norm, unit_normal, valid

    def query_batch(self, Y_local):
        g = Y_local[:, 2] - self.z_top
        raw_grad = np.tile([0.0, 0.0, 1.0], (len(Y_local), 1))
        grad_norm = np.ones(len(Y_local))
        unit_normal = np.tile([0.0, 0.0, 1.0], (len(Y_local), 1))
        valid = np.ones(len(Y_local), dtype=bool)
        return g, raw_grad, grad_norm, unit_normal, valid
