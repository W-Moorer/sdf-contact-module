import numpy as np

def analytic_torsion_torque(N, mu, R_o, R_i, omega_sign=1.0):
    if R_o == R_i:
        return 0.0
    R_spin = (2.0 / 3.0) * (R_o**3 - R_i**3) / (R_o**2 - R_i**2)
    T_z = -mu * N * R_spin * np.sign(omega_sign)
    return T_z, R_spin

def analytic_normal_force(R_o, R_i, p):
    N = p * np.pi * (R_o**2 - R_i**2)
    return N

def analytic_equivalent_radius(R_o, R_i):
    if R_o == R_i:
        return 0.0
    return (2.0 / 3.0) * (R_o**3 - R_i**3) / (R_o**2 - R_i**2)
