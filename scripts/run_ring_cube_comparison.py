"""Ring-Cube: dt=1e-6, analytic J, re-run with speedup."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import AnalyticPlaneSDF
from src.integrators.explicit_projected import ExplicitProjectedIntegrator

SC = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(SC, '..', 'models')
REF = os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision')
OUT = os.path.join(SC, '..', 'assets', 'validation')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
                    'Ring - Cube - Centric Collision_01',
                    'Ring - Cube - Centric Collision.rmd')
os.makedirs(OUT, exist_ok=True)
OUT_CSV = os.path.join(SC, '..', 'prototype', 'outputs', 'ring_cube_result.csv')

print('Loading model...')
converter = RMDToIR()
model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.quadrature_settings['damping'] = 30000.0

R_i,R_o=0.4,0.6; nr,na=8,24; dr=(R_o-R_i)/nr; dth=2*np.pi/na
Xq,wq=[],[]
for i in range(nr):
    ri=R_i+i*dr; ro=R_i+(i+1)*dr; rc=0.5*(ri+ro); a=0.5*(ro**2-ri**2)*dth
    for j in range(na):
        th=(j+0.5)*dth; Xq.append([rc*np.cos(th),rc*np.sin(th),-0.024]); wq.append(a)
qm=QuadratureMesh(np.array(Xq),np.array(wq))

eng = SDFContactEngine()
cp_obj = list(model.contacts.values())[0]
eng.register_pair(cp_obj.id, qm, AnalyticPlaneSDF(0.0),
                  aabb_b_half=[0.23,0.23,0.08], narrow_margin=0.01)
print(f'Quad: {qm.num_points} pts, k_n={cp_obj.normal_stiffness:.1e}')

print('Loading RecurDyn reference...')
rd_t, rd_z = np.loadtxt(os.path.join(REF,'Bodies_Body2_Pos_TZ.csv'),delimiter=',',skiprows=1,unpack=True)
rd_fm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_FM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_tm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_TM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_b2_vz = np.loadtxt(os.path.join(REF,'Bodies_Body2_Vel_TZ.csv'),delimiter=',',skiprows=1,unpack=True)[1]

print('Running simulation...')
state = State(model); idx4 = state.body_idx(4)
our_t, our_z, our_vz = [], [], []
intg = ExplicitProjectedIntegrator(model, contact_engine=eng)

# Phase 1: free fall
print('Phase 1: Free fall dt=0.001')
for i in range(210):
    intg.step(state, 0.001)
    our_t.append(state.t); our_z.append(state.r[idx4,2]); our_vz.append(state.v[idx4,2])

# Phase 2: contact
print('Phase 2: Contact dt=1e-6')
t0 = time.time(); sc = 0
for i in range(2000000):
    intg.step(state, 1e-6)
    our_t.append(state.t); our_z.append(state.r[idx4,2]); our_vz.append(state.v[idx4,2])
    if abs(state.v[idx4,2]) < 5e-5 and state.r[idx4,2] > 0.0238:
        sc += 1
        if sc > 2000:
            print(f'  Settled at step {i}, t={state.t:.6f}s, z={state.r[idx4,2]:.8f}')
            for _ in range(1000):
                intg.step(state, 1e-6)
                our_t.append(state.t); our_z.append(state.r[idx4,2]); our_vz.append(state.v[idx4,2])
            break
    else: sc = 0
    if i % 20000 == 0 and i > 0:
        print(f'  {i//1000:4d}k t={state.t:.4f}s z={state.r[idx4,2]:.8f} vz={state.v[idx4,2]:.4f}')
    if state.t >= rd_t[-1]: break

tot = time.time() - t0
print(f'Done: {len(our_t)} steps in {tot:.0f}s')

our_t = np.array(our_t); our_z_mm = np.array(our_z)*1000

# Save CSV
np.savetxt(OUT_CSV, np.column_stack([our_t, our_z, our_vz]),
           delimiter=',', header='time,z_m,vz_ms', comments='')
print(f'Saved: {OUT_CSV}')

# Plot
F_eq = 2.6148 * 9.80665
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax = axes[0]
ax.plot(rd_t, rd_fm, 'b-', lw=2, label='RecurDyn')
ax.plot([rd_t[0],rd_t[-1]],[F_eq,F_eq],'k:',lw=1.5,label=f'F=m*g={F_eq:.1f}N')
ax.set_title('Contact Force'); ax.legend(); ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(rd_t, rd_tm, 'r-', lw=2, label='RecurDyn')
ax.plot([rd_t[0],rd_t[-1]],[0,0],'k:',lw=1.5)
ax.set_title('Contact Torque'); ax.legend(); ax.grid(alpha=0.3)

ax = axes[2]
ax.plot(rd_t, rd_z, 'b-', lw=2, label='RecurDyn')
ax.plot(our_t, our_z_mm, 'r--', lw=2, label=f'Our (dt=1e-6)')
ax.plot([rd_t[0],rd_t[-1]],[rd_z[-1],rd_z[-1]],'k:',lw=1.5,label=f'RecurDyn eq={rd_z[-1]:.0f}mm')
our_settled = our_z[-1]*1000
ax.plot([rd_t[0],rd_t[-1]],[our_settled,our_settled],'m:',lw=1.5,label=f'Our eq={our_settled:.0f}mm')
ax.set_title('Ring Z Position'); ax.legend(); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT,'ring_cube_comparison.png'),dpi=150,bbox_inches='tight')
plt.close()
print(f'Plot: {OUT}/ring_cube_comparison.png')

print(f'\nRecurDyn Z={rd_z[-1]:.1f}mm | Our Z={our_z[-1]*1000:.1f}mm | Vz={state.v[idx4,2]:.2e}')
