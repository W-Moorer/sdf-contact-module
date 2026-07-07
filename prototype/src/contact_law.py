import numpy as np

def compute_normal_pressure(gap, k_n, epsilon=0.0):
    pressure = k_n * np.maximum(epsilon - gap, 0.0)
    return pressure

def compute_friction_traction(u_t, p_hat, mu, regularizer=0.0):
    u_t_norm = np.linalg.norm(u_t)
    if regularizer > 0.0:
        denom = np.sqrt(u_t_norm**2 + regularizer**2)
    else:
        denom = u_t_norm if u_t_norm > 1e-15 else 1.0
    f_t = -mu * p_hat * u_t / denom
    return f_t
