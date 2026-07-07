import numpy as np

def load_obj_triangles(path):
    verts = []
    tris = []
    with open(path) as f:
        for line in f:
            if line.startswith('v '):
                verts.append([float(x) for x in line.strip().split()[1:]])
            elif line.startswith('f '):
                parts = line.strip().split()
                idx = [int(p.split('/')[0]) - 1 for p in parts[1:]]
                if len(idx) == 3:
                    tris.append(idx)
    verts = np.array(verts)
    tris = np.array(tris)
    return verts, tris

def triangle_centroid_quadrature(verts, tris):
    v = verts[tris]
    X_q = v.mean(axis=1)
    edges0 = v[:, 1] - v[:, 0]
    edges1 = v[:, 2] - v[:, 0]
    normals = np.cross(edges0, edges1)
    areas = 0.5 * np.linalg.norm(normals, axis=1)
    w_q = areas.copy()
    tri_normals = normals / (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-30)
    return X_q, w_q, tri_normals, np.arange(len(tris))

def select_bottom_face(X_q, tri_normals, tol=0.01):
    mask = tri_normals[:, 2] < -(1.0 - tol)
    return mask

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
