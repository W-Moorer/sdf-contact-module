import numpy as np
from scipy.sparse import block_diag, csr_matrix


def body_inertia_world(R, inertia_body):
    return R @ inertia_body @ R.T


def build_mass_matrix(model, state):
    blocks = []
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        m = body.mass
        Iw = body_inertia_world(state.R[idx], body.inertia_body)
        block = np.zeros((6, 6))
        block[:3, :3] = m * np.eye(3)
        block[3:, 3:] = Iw
        blocks.append(block)
    if len(blocks) == 0:
        return csr_matrix((0, 0))
    return block_diag(blocks, format='csr')


def build_gyroscopic(model, state):
    nb_m = model.num_movable
    C = np.zeros(6 * nb_m)
    tmp_idx = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        omega = state.omega[idx]
        Iw = body_inertia_world(state.R[idx], body.inertia_body)
        C[6*tmp_idx + 3:6*tmp_idx + 6] = np.cross(omega, Iw @ omega)
        tmp_idx += 1
    return C
