"""Trace first contact with correct offsets, RMD params, no FD."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR; from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.integrators.implicit_be import BackwardEulerIntegrator

SC = os.path.join(os.path.dirname(__file__), '..')
RMD = os.path.join(SC, 'models', '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id

def load_q(path):
    with open(path) as f: l = f.readlines()
    v, t = [], []
    for x in l:
        if x.startswith('v '): v.append([float(y) for y in x.strip().split()[1:4]])
        elif x.startswith('f '): t.append([int(y.split('/')[0])-1 for y in x.strip().split()[1:4]])
    v = np.array(v); X, w = [], []
    for tri in t:
        vv = v[tri]; X.append(vv.mean(0)); w.append(0.5*np.linalg.norm(np.cross(vv[1]-vv[0], vv[2]-vv[0])))
    return QuadratureMesh(np.array(X), np.array(w))

qm = load_q(os.path.join(SC, 'models', 'rmd_ring.obj'))
qm.X_q += np.array([-0.099, -0.099, -0.024])
eng = SDFContactEngine(); cp_obj = list(model.contacts.values())[0]
sdf = TrilinearSDFGrid(os.path.join(SC, 'models', 'rmd_box_res128.sdf'))
eng.register_pair(cp_obj.id, qm, sdf, aabb_b_half=[0.23,0.23,0.5],
                  narrow_margin=1.0, marker_offset_b=np.array([-0.23, -0.23, -0.08]))

state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-5, contact_dt_min=1e-8,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 1e-5; intg._stable_vn_tol = 1e-3

# Run until near contact
for i in range(490):
    intg.step(state, 0.001)

# Trace contact steps
for i in range(20):
    m = eng.measure_pair(cp_obj, state, model)
    Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
    ph = list(intg._contact_states.values())[0].phase.name if intg._contact_states else '?'
    print(f'Step{i+490:3d}: t={state.t:.4f} z={state.r[idx4,2]*1000:.2f}mm vz={state.v[idx4,2]:.4f} gap={m.gap*1000:.4f}mm act={m.active_count} Fz={Q[6*1+2]:.2f}N ph={ph}')
    intg.step(state, 0.001)
    if state.r[idx4,2] < -0.1:
        print(f'  FELL THROUGH at step {i+490}')
        break
