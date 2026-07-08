"""Trace TOI with full surface quad mesh."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine, ContactMeasure
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.integrators.implicit_be import BackwardEulerIntegrator, ContactPhase
import time as _time

SC = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(SC, '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')

converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.activation_distance = 2e-3; cp.quadrature_settings['damping'] = 30000.0

def load_quad(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    verts, tris = [], []
    for line in lines:
        if line.startswith('v '):
            verts.append([float(x) for x in line.strip().split()[1:4]])
        elif line.startswith('f '):
            parts = line.strip().split()
            tris.append([int(p.split('/')[0])-1 for p in parts[1:4]])
    verts = np.array(verts)
    X_q, w_q = [], []
    for t in tris:
        v = verts[t]
        X_q.append((v[0]+v[1]+v[2])/3)
        w_q.append(0.5*np.linalg.norm(np.cross(v[1]-v[0], v[2]-v[0])))
    X_q = np.array(X_q)
    X_q[:, 2] = -X_q[:, 2]  # TOP→BOTTOM
    return QuadratureMesh(X_q, np.array(w_q))

qm = load_quad(os.path.join(MODELS, 'rmd_ring.obj'))
print(f'Quad: {qm.num_points} pts, area={np.sum(qm.w_q):.4f}m²')
print(f'  z-range: [{qm.X_q[:,2].min():.4f}, {qm.X_q[:,2].max():.4f}]')

eng = SDFContactEngine()
cp_obj = list(model.contacts.values())[0]
cube_sdf = TrilinearSDFGrid(os.path.join(MODELS, 'rmd_box_res128.sdf'))
eng.register_pair(cp_obj.id, qm, cube_sdf, aabb_b_half=[0.23,0.23,0.5], narrow_margin=0.2)

state = State(model)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4, contact_dt_min=1e-6,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 2e-3; intg._stable_vn_tol = 1e-3

idx4 = state.body_idx(4)
for i in range(139):
    intg.step(state, 0.001)

# After TOI placement
print(f'After TOI: t={state.t:.6f} z={state.r[idx4,2]:.5f} vz={state.v[idx4,2]:.4f}')
m = eng.measure_pair(cp_obj, state, model)
print(f'gap={m.gap*1e6:.0f}um act={m.active_count} valid={m.valid_count} pen={m.penetration*1e6:.0f}um')
Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
print(f'Fz={Q[6*1+2]:.1f}N')

# Take steps and trace
for i in range(10):
    intg.step(state, 0.001)
    m = eng.measure_pair(cp_obj, state, model)
    Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
    ph = list(intg._contact_states.values())[0].phase.name
    print(f'Step{140+i}: t={state.t:.4f} z={state.r[idx4,2]:.5f} vz={state.v[idx4,2]:.4f} gap={m.gap*1e6:.0f}um act={m.active_count} Fz={Q[6*1+2]:.0f}N ph={ph}')
    if state.r[idx4,2] < -0.05:
        print('  FELL THROUGH!'); break
