"""Ring-cube with exact GSurface mesh quadrature from OBJ."""
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
import time as _time

SC = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(SC, '..', 'models')
REF = os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision')
OUT = os.path.join(SC, '..', 'assets', 'validation')
os.makedirs(OUT, exist_ok=True)

# --- Load model ---
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
print('Loading model...')
converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.activation_distance = 2e-3
    cp.quadrature_settings['damping'] = 30000.0

# --- Load ring OBJ and compute triangle-centroid quadrature ---
def load_obj_quadrature(path):
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
        v0, v1, v2 = verts[t[0]], verts[t[1]], verts[t[2]]
        centroid = (v0 + v1 + v2) / 3.0
        area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
        X_q.append(centroid)
        w_q.append(area)
    X_q = np.array(X_q)
    # Transform from marker frame to body COM frame using Frame 4 offset
    # Frame 4: Body2.Torus1.BaseGSurfacePatchRefMarker, pos=[-0.099, -0.099, -0.024]
    marker_offset = np.array([-0.099, -0.099, -0.024])
    X_q = X_q + marker_offset  # R=I for this marker
    return QuadratureMesh(X_q, np.array(w_q))

ring_obj_path = os.path.join(MODELS, 'rmd_ring.obj')
print(f'Loading ring mesh: {ring_obj_path}')
qm = load_obj_quadrature(ring_obj_path)
print(f'Quad: {qm.num_points} pts, total area={np.sum(qm.w_q):.4f} m²')
print(f'  Local x: [{qm.X_q[:,0].min():.4f}, {qm.X_q[:,0].max():.4f}]')
print(f'  Local y: [{qm.X_q[:,1].min():.4f}, {qm.X_q[:,1].max():.4f}]')
print(f'  Local z: [{qm.X_q[:,2].min():.4f}, {qm.X_q[:,2].max():.4f}]')

# --- Setup contact engine ---
eng = SDFContactEngine()
cp_obj = list(model.contacts.values())[0]
cube_sdf = TrilinearSDFGrid(os.path.join(MODELS, 'rmd_box_res128.sdf'))
# Cube SDF was generated from Frame 3: Body1.Box1.BaseGSurfacePatchRefMarker
# marker_offset=[-0.23, -0.23, -0.08] — offset from cube COM to SDF coordinate frame
cube_marker_offset = np.array([-0.23, -0.23, -0.08])
eng.register_pair(cp_obj.id, qm, cube_sdf, aabb_b_half=[0.23,0.23,0.5],
                  narrow_margin=0.2, marker_offset_b=cube_marker_offset)

# --- Load reference ---
print('Loading RecurDyn reference...')
rd_t, rd_z = np.loadtxt(os.path.join(REF,'Bodies_Body2_Pos_TZ.csv'),delimiter=',',skiprows=1,unpack=True)
rd_fm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_FM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_tm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_TM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_b2_vz = np.loadtxt(os.path.join(REF,'Bodies_Body2_Vel_TZ.csv'),delimiter=',',skiprows=1,unpack=True)[1]

# --- Run ---
print('Running BE simulation...')
state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4, contact_dt_min=1e-6,
    contact_dt_impact=1e-3, contact_dt_active=1e-3, contact_dt_max=1e-3, contact_dt_normal=0.001)
intg._release_gap = 2e-3
intg._stable_vn_tol = 1e-3

t0 = _time.perf_counter()
t_arr, z_arr, vz_arr = [], [], []
def cb(state, i, total):
    t_arr.append(state.t); idx = state.body_idx(4)
    z_arr.append(state.r[idx,2].copy()); vz_arr.append(state.v[idx,2].copy())
    if len(t_arr) % 1000 == 0:
        print(f'  step{len(t_arr):4d} t={state.t:.4f} z={state.r[idx,2]:.5f} vz={state.v[idx,2]:.4f}')
intg.integrate(state, 1.0, 0.001, callback=cb)
t = _time.perf_counter() - t0
print(f'Done: {len(t_arr)} steps in {t:.1f}s ({len(t_arr)/t:.0f} steps/s)')

t_arr, z_arr, vz_arr = np.array(t_arr), np.array(z_arr), np.array(vz_arr)
z_mm = z_arr * 1000

# --- Save + Plot ---
csv_path = os.path.join(SC, '..', 'prototype', 'outputs', 'be_full_surface.csv')
np.savetxt(csv_path, np.column_stack([t_arr, z_arr, vz_arr]), delimiter=',', header='time,z_m,vz_ms', comments='')
print(f'Saved: {csv_path}')

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]; ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, z_mm, 'r--', lw=1.5, label='BE (full surface)'); ax.set_title('Ring Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[1]; ax.plot(t_arr, z_arr, 'r-', lw=1.5); ax.axhline(y=0.154, color='k', ls=':', label='cube top'); ax.set_title('BE Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[2]; ax.plot(rd_t, rd_b2_vz, 'b-', lw=1.5, label='RecurDyn'); ax.plot(t_arr, vz_arr, 'g--', lw=1.5, label='BE'); ax.set_title('Z Velocity'); ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'be_full_surface.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f'Plot saved')
print(f'Final z={z_arr[-1]*1000:.1f}mm (Ref={rd_z[-1]:.1f}mm)')
