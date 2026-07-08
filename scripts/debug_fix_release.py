"""Debug: fix release_gap and test contact stability."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.integrators.implicit_be import BackwardEulerIntegrator

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
converter = RMDToIR()
model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.activation_distance = 2e-3
    cp.quadrature_settings['damping'] = 30000.0

R_i, R_o = 0.4, 0.6; h_half = 0.024
nr, na, nz = 16, 48, 8
dr = (R_o - R_i) / nr; dth = 2 * np.pi / na; dz = 2 * h_half / nz
Xq, wq = [], []
for i in range(nr):
    ri = R_i + i*dr; ro = R_i + (i+1)*dr; rc = 0.5*(ri+ro)
    a = 0.5*(ro**2 - ri**2)*dth
    for j in range(na):
        th = (j+0.5)*dth; Xq.append([rc*np.cos(th), rc*np.sin(th), -h_half]); wq.append(a)
R_in = R_i; A_cyl = dth * R_in * dz
for k in range(nz):
    z_pt = -h_half + (k+0.5)*dz
    for j in range(na):
        th = (j+0.5)*dth; Xq.append([R_in*np.cos(th), R_in*np.sin(th), z_pt]); wq.append(A_cyl)
R_out = R_o
for k in range(nz):
    z_pt = -h_half + (k+0.5)*dz
    for j in range(na):
        th = (j+0.5)*dth; Xq.append([R_out*np.cos(th), R_out*np.sin(th), z_pt]); wq.append(A_cyl)
qm = QuadratureMesh(np.array(Xq), np.array(wq))
print(f'Quad: {qm.num_points} pts')

eng = SDFContactEngine()
cp_obj = list(model.contacts.values())[0]
cube_sdf = TrilinearSDFGrid(os.path.join(MODELS, 'rmd_box_res128.sdf'))
eng.register_pair(cp_obj.id, qm, cube_sdf, aabb_b_half=[0.23,0.23,0.5], narrow_margin=0.2)

state = State(model)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4, contact_dt_min=1e-6,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
# Fix: release_gap should equal activation_distance so contact stays active within zone
intg._release_gap = 2e-3
intg._stable_vn_tol = 1e-3  # Relax velocity tolerance

idx4 = state.body_idx(4)
t0 = time.perf_counter() if 'time' in dir() else __import__('time').perf_counter
import time as _time
t0 = _time.perf_counter()
for i in range(160):
    intg.step(state, 0.001)
    if i >= 137:
        m = eng.measure_pair(cp_obj, state, model)
        Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
        ph = list(intg._contact_states.values())[0].phase.name
        print(f'Step{i:3d}: t={state.t:.4f} z={state.r[idx4,2]:.5f} vz={state.v[idx4,2]:.4f} gap={m.gap*1e6:.0f}um act={m.active_count} Fz={Q[6*1+2]:.0f}N ph={ph}')
    if state.r[idx4,2] < -0.05 and i > 145:
        print('  FELL THROUGH!'); break
t = _time.perf_counter() - t0
print(f'Time: {t:.1f}s')
