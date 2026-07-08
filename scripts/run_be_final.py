"""Final run: ring-cube with real SDF, enhanced quad mesh, corrected release_gap."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.integrators.implicit_be import BackwardEulerIntegrator

SC = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(SC, '..', 'models')
REF = os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision')
OUT = os.path.join(SC, '..', 'assets', 'validation')
os.makedirs(OUT, exist_ok=True)

RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
print('Loading model...')
converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.activation_distance = 2e-3
    cp.quadrature_settings['damping'] = 30000.0

R_i, R_o = 0.4, 0.6; h_half = 0.024; nr, na, nz = 16, 48, 8
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

print('Loading RecurDyn reference...')
rd_t, rd_z = np.loadtxt(os.path.join(REF,'Bodies_Body2_Pos_TZ.csv'),delimiter=',',skiprows=1,unpack=True)
rd_fm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_FM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_tm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_TM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_b2_vz = np.loadtxt(os.path.join(REF,'Bodies_Body2_Vel_TZ.csv'),delimiter=',',skiprows=1,unpack=True)[1]

print('Running BE simulation...')
state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4, contact_dt_min=1e-6,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 2e-3  # match activation_distance
intg._stable_vn_tol = 1e-3

import time as _time
t0 = _time.perf_counter()
t_arr, z_arr, vz_arr = [], [], []
def cb(state,i,total):
    t_arr.append(state.t); idx=state.body_idx(4)
    z_arr.append(state.r[idx,2].copy()); vz_arr.append(state.v[idx,2].copy())
    if len(t_arr) % 1000 == 0:
        print(f'  step{len(t_arr):4d} t={state.t:.4f} z={state.r[idx,2]:.5f} vz={state.v[idx,2]:.4f}')
intg.integrate(state, 0.5, 0.001, callback=cb)
t = _time.perf_counter() - t0
print(f'Done: {len(t_arr)} steps in {t:.1f}s ({len(t_arr)/t:.0f} steps/s)')

t_arr, z_arr, vz_arr = np.array(t_arr), np.array(z_arr), np.array(vz_arr)
z_mm = z_arr * 1000

# Save CSV
csv_path = os.path.join(SC, '..', 'prototype', 'outputs', 'be_final.csv')
np.savetxt(csv_path, np.column_stack([t_arr, z_arr, vz_arr]), delimiter=',', header='time,z_m,vz_ms', comments='')
print(f'Saved: {csv_path}')

# Plot
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]; ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, z_mm, 'r--', lw=1.5, label='BE (SDF)'); ax.set_title('Ring Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[1]; ax.plot(t_arr, z_arr, 'r-', lw=1.5); ax.axhline(y=0.154, color='k', ls=':', label='contact'); ax.set_title('BE Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[2]; ax.plot(rd_t, rd_b2_vz, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, vz_arr, 'g--', lw=1.5, label='BE'); ax.set_title('Z Velocity'); ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'be_final.png'), dpi=150, bbox_inches='tight')
plt.close()

print(f'Final z={z_arr[-1]*1000:.1f}mm (Ref={rd_z[-1]:.1f}mm)')
print(f'Min z={z_arr.min()*1000:.1f}mm Max z={z_arr.max()*1000:.1f}mm')
