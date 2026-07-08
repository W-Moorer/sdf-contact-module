"""Ring-Cube collision with adaptive-contact BE integrator (prototype_adaptive)."""
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
from src.sdf_grid import AnalyticPlaneSDF
from src.integrators.implicit_be import BackwardEulerIntegrator

SC = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(SC, '..', 'models')
REF = os.path.join(SC, '..', 'reference_csv', 'ring_cube_collision')
OUT = os.path.join(SC, '..', 'assets', 'validation')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
                    'Ring - Cube - Centric Collision_01',
                    'Ring - Cube - Centric Collision.rmd')
os.makedirs(OUT, exist_ok=True)

print('Loading model...')
converter = RMDToIR()
model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.quadrature_settings['damping'] = 30000.0
    cp.activation_distance = 2e-4

R_i,R_o=0.4,0.6; nr,na=8,24; dr=(R_o-R_i)/nr; dth=2*np.pi/na
Xq,wq=[],[]
for i in range(nr):
    ri=R_i+i*dr;ro=R_i+(i+1)*dr;rc=0.5*(ri+ro);a=0.5*(ro**2-ri**2)*dth
    for j in range(na):
        th=(j+0.5)*dth;Xq.append([rc*np.cos(th),rc*np.sin(th),-0.024]);wq.append(a)
qm=QuadratureMesh(np.array(Xq),np.array(wq))
eng=SDFContactEngine()
eng._narrow_margin = 0.2  # Ensure narrow filter doesn't block contact
cp_obj=list(model.contacts.values())[0]
eng.register_pair(cp_obj.id, qm, AnalyticPlaneSDF(0.13),
                  aabb_b_half=[0.23,0.23,0.5], narrow_margin=0.2)

print('Loading RecurDyn reference...')
rd_t, rd_z = np.loadtxt(os.path.join(REF,'Bodies_Body2_Pos_TZ.csv'),delimiter=',',skiprows=1,unpack=True)
rd_fm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_FM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_tm = np.loadtxt(os.path.join(REF,'Contact_Geo_Contact_GeoSurContact1_TM_GeoContact.csv'),delimiter=',',skiprows=1,unpack=True)[1]
rd_b2_vz = np.loadtxt(os.path.join(REF,'Bodies_Body2_Vel_TZ.csv'),delimiter=',',skiprows=1,unpack=True)[1]

print('Running BE simulation with adaptive contact stepping...')
state = State(model); idx4 = state.body_idx(4)
intg = BackwardEulerIntegrator(
    model, contact_engine=eng, max_iter=50, tol=1e-8,
    contact_penetration_tol=1e-4,
    contact_dt_min=1e-6,
    contact_dt_impact=1e-3,
    contact_dt_active=1e-3,
    contact_dt_max=1e-3,
    contact_dt_normal=0.001,
)

t_arr, z_arr, vz_arr, nsteps = [], [], [], 0
print_record = [0.0]

def callback(state, i, total):
    t_arr.append(state.t); z_arr.append(state.r[idx4,2].copy()); vz_arr.append(state.v[idx4,2].copy())
    if state.t - print_record[0] >= 0.1:
        print(f'  t={state.t:.4f}s z={state.r[idx4,2]:.5f}m vz={state.v[idx4,2]:.4f} '
              f'phases={[st.phase.name for st in intg._contact_states.values()]}')
        print_record[0] = state.t

intg.integrate(state, 1.0, 0.001, callback=callback)
nsteps = len(t_arr)

print(f'Completed: {nsteps} steps, t={state.t:.4f}s, z={state.r[idx4,2]:.4f}m')
t_arr, z_arr, vz_arr = np.array(t_arr), np.array(z_arr), np.array(vz_arr)
z_mm = z_arr * 1000

# Save CSV
csv_path = os.path.join(SC, '..', 'prototype', 'outputs', 'be_ring_cube_result.csv')
np.savetxt(csv_path, np.column_stack([t_arr, z_arr, vz_arr]),
           delimiter=',', header='time,z_m,vz_ms', comments='')
print(f'Saved: {csv_path}')

# Plot comparison
fig, axes = plt.subplots(2, 3, figsize=(18, 8))
ax = axes[0,0]; ax.plot(rd_t, rd_fm, 'b-', lw=1.5, label='RecurDyn'); ax.set_title('Contact Force (RecurDyn)'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[0,1]; ax.plot(rd_t, rd_tm, 'r-', lw=1.5, label='RecurDyn'); ax.set_title('Contact Torque (RecurDyn)'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[0,2]; ax.plot(rd_t, rd_z, 'b-', lw=1.5, label='RecurDyn')
ax.plot(t_arr, z_mm, 'r--', lw=1.5, label=f'BE adaptive')
ax.set_title('Ring Z Position'); ax.legend(); ax.grid(alpha=0.3)

ax = axes[1,0]; ax.plot(t_arr, z_arr, 'r-', lw=1.5); ax.axhline(y=0.154, color='k', ls=':', label='contact z=0.154')
ax.set_title('BE Z Position'); ax.legend(); ax.grid(alpha=0.3)
ax = axes[1,1]; ax.plot(t_arr, vz_arr, 'g-', lw=1.5); ax.axhline(y=0, color='k', ls=':')
ax.set_title('BE Z Velocity'); ax.grid(alpha=0.3)
ax = axes[1,2]; ax.plot(rd_t, rd_b2_vz, 'b-', lw=1.5, label='RecurDyn')
ax.plot(t_arr, vz_arr, 'g--', lw=1.5, label='BE'); ax.set_title('Z Velocity Comparison'); ax.legend(); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUT,'be_ring_cube_comparison.png'),dpi=150,bbox_inches='tight')
plt.close()
print(f'Plot: {OUT}/be_ring_cube_comparison.png')

if len(z_arr) > 0:
    print(f'Final: z={z_arr[-1]*1000:.1f}mm (Ref: {rd_z[-1]:.1f}mm)')
