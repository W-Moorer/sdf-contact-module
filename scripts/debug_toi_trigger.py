"""Debug TOI trigger with real SDF + enhanced quadrature."""
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
    ri = R_i + i * dr; ro = R_i + (i + 1) * dr; rc = 0.5 * (ri + ro)
    a = 0.5 * (ro**2 - ri**2) * dth
    for j in range(na):
        th = (j + 0.5) * dth
        Xq.append([rc * np.cos(th), rc * np.sin(th), -h_half]); wq.append(a)
R_in = R_i; A_cyl = dth * R_in * dz
for k in range(nz):
    z_pt = -h_half + (k + 0.5) * dz
    for j in range(na):
        th = (j + 0.5) * dth; Xq.append([R_in * np.cos(th), R_in * np.sin(th), z_pt]); wq.append(A_cyl)
R_out = R_o
for k in range(nz):
    z_pt = -h_half + (k + 0.5) * dz
    for j in range(na):
        th = (j + 0.5) * dth; Xq.append([R_out * np.cos(th), R_out * np.sin(th), z_pt]); wq.append(A_cyl)

qm = QuadratureMesh(np.array(Xq), np.array(wq))
eng = SDFContactEngine()
cp_obj = list(model.contacts.values())[0]
cube_sdf = TrilinearSDFGrid(os.path.join(MODELS, 'rmd_box_res128.sdf'))
eng.register_pair(cp_obj.id, qm, cube_sdf, aabb_b_half=[0.23, 0.23, 0.5], narrow_margin=0.2)

state = State(model)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4, contact_dt_min=1e-6,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)

idx4 = state.body_idx(4)
for i in range(138):
    intg.step(state, 0.001)

print(f'Step 138: t={state.t:.4f} z={state.r[idx4,2]:.6f} vz={state.v[idx4,2]:.4f}')
m = eng.measure_pair(cp_obj, state, model)
print(f'  gap={m.gap*1e6:.0f}um active={m.active_count} valid={m.valid_count}')

# Manual trial
print('\nTrial step at dt=0.001 (no contact forces)...')
trial = state.copy()
intg._step_be(trial, 0.001, record=False, active_contact_ids=set())
print(f'  Trial: t={trial.t:.4f} z={trial.r[idx4,2]:.6f} vz={trial.v[idx4,2]:.4f}')
m_t = eng.measure_pair(cp_obj, trial, model)
print(f'  gap={m_t.gap*1e6:.0f}um active={m_t.active_count}')

# Check if crossing detected
from src.integrators.implicit_be import ContactPhase
inactive_ids = [cid for cid, st in intg._contact_states.items() if st.phase == ContactPhase.INACTIVE]
active_ids = intg._active_contact_ids()
start_measures = intg._measure_contacts(state, cp_ids=inactive_ids)
end_measures = intg._measure_contacts(trial, cp_ids=inactive_ids)

for cid in inactive_ids:
    m0 = start_measures.get(cid); m1 = end_measures.get(cid)
    if m0 and m1:
        crosses = (m0.gap > 0.0 and m1.gap <= 0.0)
        print(f'\nPair {cid}: start_gap={m0.gap*1e6:.0f}um trial_gap={m1.gap*1e6:.0f}um crosses={crosses}')

# One more step
h_done = intg.step(state, 0.001)
print(f'\nStep 139: h_done={h_done:.6f} t={state.t:.6f} z={state.r[idx4,2]:.6f}')
print(f'Phases: {[st.phase.name for st in intg._contact_states.values()]}')
m = eng.measure_pair(cp_obj, state, model)
print(f'gap={m.gap*1e6:.0f}um active={m.active_count}')

# Step 140
h_done = intg.step(state, 0.001)
print(f'\nStep 140: h_done={h_done:.6f} t={state.t:.6f} z={state.r[idx4,2]:.6f} vz={state.v[idx4,2]:.4f}')
print(f'Phases: {[st.phase.name for st in intg._contact_states.values()]}')
Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
print(f'Contact Fz={Q[6*1+2]:.1f}N')
m = eng.measure_pair(cp_obj, state, model)
print(f'gap={m.gap*1e6:.0f}um active={m.active_count}')

# Step 141
h_done = intg.step(state, 0.001)
print(f'\nStep 141: h_done={h_done:.6f} t={state.t:.6f} z={state.r[idx4,2]:.6f} vz={state.v[idx4,2]:.4f}')
print(f'Phases: {[st.phase.name for st in intg._contact_states.values()]}')
Q = eng.compute_Q_contact(model, state, active_ids={cp_obj.id})
print(f'Contact Fz={Q[6*1+2]:.1f}N')
m = eng.measure_pair(cp_obj, state, model)
print(f'gap={m.gap*1e6:.0f}um active={m.active_count} valid={m.valid_count}')
