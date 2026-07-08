"""Ring-cube: exact RMD parameters, no tuning. KORDER=2, K=1e8, C=10, BPEN=1e-5."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np; import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import time as _time
from src.importers import RMDToIR; from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.integrators.implicit_be import BackwardEulerIntegrator

SC = os.path.dirname(os.path.abspath(__file__))
RMD = os.path.join(SC, '..', 'models', '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    # Exact RMD parameters (importer already sets them):
    # normal_stiffness=1e8, k_order=2, activation_distance=1e-5, damping=10

def load_quad(path):
    with open(path) as f: lines = f.readlines()
    verts, tris = [], []
    for line in lines:
        if line.startswith('v '): verts.append([float(x) for x in line.strip().split()[1:4]])
        elif line.startswith('f '): tris.append([int(p.split('/')[0])-1 for p in line.strip().split()[1:4]])
    verts = np.array(verts); X_q, w_q = [], []
    for t in tris:
        v = verts[t]; X_q.append(v.mean(0)); w_q.append(0.5 * np.linalg.norm(np.cross(v[1]-v[0], v[2]-v[0])))
    return QuadratureMesh(np.array(X_q), np.array(w_q))

qm = load_quad(os.path.join(SC, '..', 'models', 'rmd_ring.obj'))
eng = SDFContactEngine(); cp_obj = list(model.contacts.values())[0]
sdf = TrilinearSDFGrid(os.path.join(SC, '..', 'models', 'rmd_box_res128.sdf'))
eng.register_pair(cp_obj.id, qm, sdf, aabb_b_half=[0.23,0.23,0.5], narrow_margin=0.2)

state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-5, contact_dt_min=1e-8,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 1e-5; intg._stable_vn_tol = 1e-3

print('Running...')
t0 = _time.perf_counter(); t_arr, z_arr, vz_arr = [], [], []
def cb(s, i, t):
    t_arr.append(s.t); z_arr.append(s.r[idx4,2].copy()); vz_arr.append(s.v[idx4,2].copy())
    if len(t_arr) % 500 == 0: print(f'  step{len(t_arr):4d} t={s.t:.4f} z={s.r[idx4,2]*1000:.1f}mm vz={s.v[idx4,2]:.4f}')
intg.integrate(state, 0.5, 0.001, callback=cb)
t = _time.perf_counter() - t0
t_arr, z_arr, vz_arr = np.array(t_arr), np.array(z_arr), np.array(vz_arr)
print(f'{len(t_arr)} steps in {t:.1f}s')

# Plot
rd_t, rd_z = np.loadtxt(os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision', 'Bodies_Body2_Pos_TZ.csv'), delimiter=',', skiprows=1, unpack=True)
rd_vz = np.loadtxt(os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision', 'Bodies_Body2_Vel_TZ.csv'), delimiter=',', skiprows=1, unpack=True)[1]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]; ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, z_arr*1000, 'r--', lw=1.2, label='BE (RMD params)'); ax.axhline(y=130, color='k', ls=':', lw=0.8); ax.set_xlabel('t(s)'); ax.set_ylabel('Z(mm)'); ax.set_title('Ring Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[1]; ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, z_arr*1000, 'r--', lw=1.2, label='BE'); ax.axhline(y=130, color='k', ls=':'); ax.set_xlim(0.12, 0.35); ax.set_ylim(60, 170); ax.set_xlabel('t(s)'); ax.set_title('Contact Phase'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[2]; ax.plot(rd_t, rd_vz, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, vz_arr, 'g--', lw=1.2, label='BE'); ax.axhline(y=0, color='k', ls=':'); ax.set_xlabel('t(s)'); ax.set_ylabel('Vz(m/s)'); ax.set_title('Z Velocity'); ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(SC, '..', 'assets', 'validation', 'be_rmd_exact.png'), dpi=150); plt.close()
print(f'Plot saved. Final z={z_arr[-1]*1000:.1f}mm (Ref={rd_z[-1]:.1f}mm)')
