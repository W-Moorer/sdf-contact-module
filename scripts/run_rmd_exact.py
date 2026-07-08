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
    cp.action_marker_id, cp.base_marker_id = cp.base_marker_id, cp.action_marker_id
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
# Quad points stay in GSurface RM marker frame (action_marker_id from RMD)
eng = SDFContactEngine(); cp_obj = list(model.contacts.values())[0]
sdf = TrilinearSDFGrid(os.path.join(SC, '..', 'models', 'rmd_box_res128.sdf'))
# SDF query in base marker frame via action_frame.relative_to(base_frame, state)
eng.register_pair(cp_obj.id, qm, sdf, aabb_b_half=[0.23,0.23,0.5],
                  narrow_margin=1.0)

state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-5, contact_dt_min=1e-8,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 1e-5; intg._stable_vn_tol = 1e-3

print('Running...')
t0 = _time.perf_counter(); t_arr, z_arr, vz_arr, fz_arr = [], [], [], []
def cb(s, i, t):
    t_arr.append(s.t); z_arr.append(s.r[idx4,2].copy()); vz_arr.append(s.v[idx4,2].copy())
    Q = eng.compute_Q_contact(model, s, active_ids={cp_obj.id})
    fz_arr.append(Q[6*1 + 2])
    if len(t_arr) % 500 == 0: print(f'  step{len(t_arr):4d} t={s.t:.4f} z={s.r[idx4,2]*1000:.1f}mm vz={s.v[idx4,2]:.4f}')
intg.integrate(state, 1.0, 0.001, callback=cb)
t = _time.perf_counter() - t0
t_arr, z_arr, vz_arr, fz_arr = np.array(t_arr), np.array(z_arr), np.array(vz_arr), np.array(fz_arr)
print(f'{len(t_arr)} steps in {t:.1f}s')
z_mm = z_arr * 1000

# Plot
rd_t, rd_z = np.loadtxt(os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision', 'Bodies_Body2_Pos_TZ.csv'), delimiter=',', skiprows=1, unpack=True)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Left: Z position full time
ax = axes[0]
ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn')
ax.plot(t_arr, z_mm, 'r--', lw=1.2, label='BE (RMD params)')
ax.axhline(y=65, color='b', ls=':', lw=1.0, alpha=0.5, label='Eq (RD 65mm)')
ax.axhline(y=64.8, color='r', ls=':', lw=1.0, alpha=0.5, label='Eq (BE 64.8mm)')
ax.set_xlabel('t(s)'); ax.set_ylabel('Z (mm)'); ax.set_title('Ring Z Position')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Middle: Contact phase zoom
ax = axes[1]
ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn')
ax.plot(t_arr, z_mm, 'r--', lw=1.2, label='BE')
ax.axhline(y=65, color='b', ls=':', lw=1.0, alpha=0.5)
ax.axhline(y=64.8, color='r', ls=':', lw=1.0, alpha=0.5)
ax.set_xlim(0.12, 0.35); ax.set_ylim(60, 170)
ax.set_xlabel('t(s)'); ax.set_ylabel('Z (mm)'); ax.set_title('Contact Phase (zoom)')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Right: Contact force Fz
ax = axes[2]
ax.plot(t_arr, fz_arr, 'r-', lw=1.2, label='BE Contact Fz')
ax.axhline(y=25.6, color='k', ls=':', lw=1.0, alpha=0.7, label='Gravity mg=25.6N')
ax.set_xlabel('t(s)'); ax.set_ylabel('Contact Force Fz (N)'); ax.set_title('Contact Force')
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_ylim(0, 500)

plt.tight_layout()
plt.savefig(os.path.join(SC, '..', 'assets', 'validation', 'be_rmd_exact.png'), dpi=150); plt.close()
print(f'Plot saved. Final z={z_mm[-1]:.1f}mm (Ref={rd_z[-1]:.1f}mm), final Fz={fz_arr[-1]:.1f}N')
