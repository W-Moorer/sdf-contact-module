import numpy as np

def annular_quadrature(R_i, R_o, n_radial=20, n_azimuthal=64):
    dr = (R_o - R_i) / n_radial
    dtheta = 2.0 * np.pi / n_azimuthal

    X_q = []
    w_q = []
    tri_id = []

    idx = 0
    for i in range(n_radial):
        r_inner = R_i + i * dr
        r_outer = R_i + (i + 1) * dr
        r_center = 0.5 * (r_inner + r_outer)
        area = 0.5 * (r_outer**2 - r_inner**2) * dtheta

        for j in range(n_azimuthal):
            theta = (j + 0.5) * dtheta
            x = r_center * np.cos(theta)
            y = r_center * np.sin(theta)
            X_q.append([x, y, 0.0])
            w_q.append(area)
            tri_id.append(idx)
            idx += 1

    return np.array(X_q), np.array(w_q), np.array(tri_id)

def bottom_face_selector(X_q, tri_normals, face_centers):
    mask = np.abs(tri_normals[:, 2] + 1.0) < 1e-6
    return mask
