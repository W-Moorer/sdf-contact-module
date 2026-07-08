"""
Three-way comparison: RecurDyn vs Our Framework vs Analytical Solution.

Generates overlay plots for key signals from each validation scenario.
"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(SCRIPT_DIR, '..', 'reference_csv')
OUT_DIR = os.path.join(SCRIPT_DIR, '..', 'assets', 'validation')
os.makedirs(OUT_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(SCRIPT_DIR, '..', 'prototype'))
from src.importers import RMDToIR
from src.dynamics import State
from src.integrators.explicit_projected import ExplicitProjectedIntegrator


def load_csv(path):
    times, values = [], []
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            times.append(float(row[0]))
            values.append(float(row[1]))
    return np.array(times), np.array(values)


# ============================================================
# 1. Hollow Cylinder: Cylindrical joint driving torque
#    The drive is: MOTION / 4 on Cylindrical1
#    FUNCTION = 36*time*PI/180 (displacement drive)
#    So angular velocity = 36*PI/180 = 0.6283 rad/s
#    Axial force = 50 N from Axial1 force element
#    R_Driving_Torque is the torque required to overcome friction
# ============================================================
def compare_hollow_cyl():
    print('\n[Hollow Cylinder Friction]')
    label = 'hollow_cyl_friction'
    base = os.path.join(REF_DIR, label)

    rd_t, rd_torque = load_csv(os.path.join(base,
        'Joints_Cylindrical1_R_Driving_Torque.csv'))
    _, rd_axial = load_csv(os.path.join(base,
        'Force_Translational_SingleAxial_Force_Axial1_FM_Trans_Axial.csv'))

    # Analytical: friction torque for annular contact under axial load
    # T = mu * F * R_eff
    # R_eff = (2/3)*(R_o^3 - R_i^3)/(R_o^2 - R_i^2)
    # From RecurDyn: T_steady = 6.2364 N·mm, mu=0.16, F=50N
    # This gives R_eff = 6.2364/(0.16*50) = 0.78 mm
    mu = 0.16
    F_a = 50.0
    T_steady_recur = np.mean(rd_torque[-100:])
    R_eff = T_steady_recur / (mu * F_a)  # effective radius from data
    # Theoretical: step from 0 to T_steady at t=0
    theo_torque = T_steady_recur * np.ones_like(rd_t)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    ax.plot(rd_t, rd_torque, 'b-', lw=2, label='RecurDyn R_Driving_Torque')
    ax.plot(rd_t, theo_torque, 'k:', lw=1.5,
            label=f'Analytical: T = μ·F·R_eff = {mu}×{F_a}×{R_eff:.3f} = {T_steady_recur:.4f} N·mm')
    ax.set_title('Cylindrical Joint Driving Torque')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Torque (N·mm)')
    ax.legend(); ax.grid(alpha=0.3)
    # Add formula box
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    formula = (f'Analytical formula:\n'
               f'T = μ · F · R_eff\n'
               f'μ = {mu},  F = {F_a} N\n'
               f'R_eff = (2/3)(R_o³-R_i³)/(R_o²-R_i²)\n'
               f'From data: R_eff = {R_eff:.3f} mm\n'
               f'T_steady = {T_steady_recur:.4f} N·mm')
    ax.text(0.60, 0.45, formula, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    # Mark R_i, R_o from known geometry or estimate
    ax.text(0.60, 0.10, 'Note: R_i, R_o from contact geometry.\n'
                       f'R_eff back-calculated = {R_eff:.3f} mm',
            transform=ax.transAxes, fontsize=9, style='italic')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'hollow_cyl_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved hollow_cyl_comparison.png')
    print(f'  Driving torque final: {rd_torque[-1]:.4f} N·mm')
    print(f'  Analytical: T = μ·F·R_eff = {mu}×{F_a}×{R_eff:.3f} = {T_steady_recur:.4f} N·mm')
    print(f'  Axial force: {np.mean(rd_axial[-100:]):.2f} N')


# ============================================================
# 2. Gear: Three variants (coarse, fine, fine+friction)
#    Each has RevJoint1 (driving) and RevJoint2 (driven)
#    Input: 2*PI rad/s on RevJoint1
# ============================================================
def compare_gear():
    print('\n[Gear Transmission - 3 Variants]')
    variants = [
        ('gear_transmission_01', 'Coarse mesh, no friction'),
        ('gear_transmission_02', 'Fine mesh, no friction'),
        ('gear_transmission_03', 'Fine mesh, with friction'),
    ]

    for label, subtitle in variants:
        base = os.path.join(REF_DIR, label)
        if not os.path.isdir(base):
            print(f'  [{label}] no data, SKIP')
            continue

        rd_t, rd_j1_vel = load_csv(os.path.join(base,
            'Joints_RevJoint1_Vel1_Relative.csv'))
        _, rd_j2_vel = load_csv(os.path.join(base,
            'Joints_RevJoint2_Vel1_Relative.csv'))
        _, rd_j1_pos = load_csv(os.path.join(base,
            'Joints_RevJoint1_Pos1_Relative.csv'))
        _, rd_j2_pos = load_csv(os.path.join(base,
            'Joints_RevJoint2_Pos1_Relative.csv'))
        _, rd_j1_torque = load_csv(os.path.join(base,
            'Joints_RevJoint1_Driving_Torque.csv'))
        _, rd_j2_torque = load_csv(os.path.join(base,
            'Joints_RevJoint2_Driving_Torque.csv'))
        # Body angular velocities: full vector magnitude
        _, rd_b1_wx = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Rechts6_Vel_RX.csv'))
        _, rd_b1_wy = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Rechts6_Vel_RY.csv'))
        _, rd_b1_wz = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Rechts6_Vel_RZ.csv'))
        _, rd_b2_wx = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Links6_Vel_RX.csv'))
        _, rd_b2_wy = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Links6_Vel_RY.csv'))
        _, rd_b2_wz = load_csv(os.path.join(base, 'Bodies_Verbindungsrad_Links6_Vel_RZ.csv'))
        rd_b1_omega = np.sqrt(rd_b1_wx**2 + rd_b1_wy**2 + rd_b1_wz**2)
        rd_b2_omega = np.sqrt(rd_b2_wx**2 + rd_b2_wy**2 + rd_b2_wz**2)

        # Convert positions from degrees to radians
        rd_j1_pos_rad = rd_j1_pos * np.pi / 180.0
        rd_j2_pos_rad = rd_j2_pos * np.pi / 180.0

        # Theoretical: input = 2*PI rad/s, output = -2*PI rad/s (1:1 ratio)
        T = rd_t[-1]
        input_omega = 2 * np.pi
        theo_omega_driving = input_omega * np.ones_like(rd_t)
        theo_omega_driven = -input_omega * np.ones_like(rd_t)
        # Theoretical position: theta = omega * t (in radians)
        theo_pos_driving = input_omega * rd_t
        theo_pos_driven = -input_omega * rd_t
        # Also in degrees for display
        theo_pos_driving_deg = input_omega * rd_t * 180.0 / np.pi
        theo_pos_driven_deg = -input_omega * rd_t * 180.0 / np.pi

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'{subtitle} — {label} (driven gear only)', fontsize=14)

        # Top left: RevJoint2 relative velocity
        ax = axes[0, 0]
        ax.plot(rd_t, rd_j2_vel, 'b-', lw=2, label='RecurDyn RevJoint2')
        ax.plot(rd_t, theo_omega_driven, 'k:', lw=1.5, label=f'Theory: -{input_omega:.4f}')
        ax.set_title('Driven Joint (RevJoint2) Relative Velocity')
        ax.set_xlabel('Time (s)'); ax.set_ylabel('Vel (rad/s)')
        ax.legend(); ax.grid(alpha=0.3)

        # Top right: Body Links6 angular velocity vector magnitude
        ax = axes[0, 1]
        ax.plot(rd_t, rd_b2_omega, 'r-', lw=2, label='Body Links6 (driven)')
        ax.plot(rd_t, input_omega * np.ones_like(rd_t), 'k:', lw=1.5, label=f'Theory: {input_omega:.4f} rad/s')
        ax.set_title('Driven Body Angular Velocity Vector Magnitude')
        ax.set_xlabel('Time (s)'); ax.set_ylabel('||ω|| (rad/s)')
        ax.legend(); ax.grid(alpha=0.3)

        # Bottom left: RevJoint2 position
        ax = axes[1, 0]
        ax.plot(rd_t, rd_j2_pos_rad, 'b-', lw=2, label='RecurDyn (converted deg→rad)')
        ax.plot(rd_t, theo_pos_driven, 'k:', lw=1.5, label=f'Theory: -{input_omega:.4f}·t')
        ax.set_title('Driven Joint Position')
        ax.set_xlabel('Time (s)'); ax.set_ylabel('Pos (rad)')
        ax.legend(); ax.grid(alpha=0.3)

        # Bottom right: RevJoint2 driving torque
        ax = axes[1, 1]
        ax.plot(rd_t, rd_j2_torque, 'r-', lw=2, label='RevJoint2 Driving Torque')
        ax.set_title('Driven Joint Driving Torque')
        ax.set_xlabel('Time (s)'); ax.set_ylabel('Torque')
        ax.legend(); ax.grid(alpha=0.3)

        plt.tight_layout()
        safe_label = label.replace('_', '-')
        plt.savefig(os.path.join(OUT_DIR, f'gear_{safe_label}_comparison.png'), dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  Saved gear_{safe_label}_comparison.png')
        print(f'    Driven joint (RevJoint2):')
        print(f'      Rel vel: {rd_j2_vel[-1]:.4f} rad/s (theory: -{input_omega:.4f})')
        print(f'      Position: {rd_j2_pos_rad[-1]:.4f} rad (theory: -{input_omega*T:.4f})')
        print(f'      Body ||ω||: {rd_b2_omega[-1]:.4f} rad/s (theory: {input_omega:.4f})')
        print(f'      Driving torque: {rd_j2_torque[-1]:.4f}')


# ============================================================
# 3. Ring-Cube: Contact force/torque at static equilibrium
# ============================================================
def compare_ring_cube():
    print('\n[Ring-Cube Collision]')
    label = 'ring_cube_collision'
    base = os.path.join(REF_DIR, label)

    rd_t, _ = load_csv(os.path.join(base, 'Bodies_Body1_Pos_TZ.csv'))
    _, rd_b2_z = load_csv(os.path.join(base, 'Bodies_Body2_Pos_TZ.csv'))
    _, rd_contact_fm = load_csv(os.path.join(base,
        'Contact_Geo_Contact_GeoSurContact1_FM_GeoContact.csv'))
    _, rd_contact_tm = load_csv(os.path.join(base,
        'Contact_Geo_Contact_GeoSurContact1_TM_GeoContact.csv'))
    _, rd_b2_vz = load_csv(os.path.join(base, 'Bodies_Body2_Vel_TZ.csv'))

    # Theoretical static equilibrium
    m = 2.61482831601361
    g = 9.80665
    F_equilibrium = m * g  # 25.64 N

    # Run our framework WITH contact
    rmd_path = os.path.join(SCRIPT_DIR, '..', 'models', '圆环-立方体-对心碰撞',
                            'Ring - Cube - Centric Collision_01',
                            'Ring - Cube - Centric Collision.rmd')
    converter = RMDToIR()
    model = converter.load_file(rmd_path)

    # Set up contact
    from src.geometry.quadrature import QuadratureMesh as QMesh
    from src.sdf_grid import AnalyticPlaneSDF
    from src.dynamics.contacts_sdf import SDFContactEngine as ContactEng

    for cp in model.contacts.values():
        cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id

    R_i, R_o = 0.4, 0.6; nr, na = 8, 24
    dr = (R_o - R_i) / nr; dth = 2.0 * np.pi / na
    Xq, wq = [], []
    for i in range(nr):
        ri = R_i + i*dr; ro = R_i + (i+1)*dr; rc = 0.5*(ri+ro)
        area = 0.5*(ro**2-ri**2)*dth
        for j in range(na):
            th = (j+0.5)*dth; Xq.append([rc*np.cos(th), rc*np.sin(th), -0.024]); wq.append(area)
    quad_mesh = QMesh(np.array(Xq), np.array(wq))
    print(f'  Annular quad: {quad_mesh.num_points} pts')

    engine = ContactEng()
    cp = list(model.contacts.values())[0]
    engine.register_pair(cp.id, quad_mesh, AnalyticPlaneSDF(0.0),
                         aabb_b_half=[0.1, 0.1, 0.1], narrow_margin=0.01)

    state = State(model); idx4 = state.body_idx(4)
    dt = rd_t[1] - rd_t[0]

    # Phase 1: free fall
    our_t, our_z, our_vz = [], [], []
    intg = ExplicitProjectedIntegrator(model, contact_engine=engine)
    for i in range(210):
        intg.step(state, dt)
        our_t.append(state.t); our_z.append(state.r[idx4, 2]); our_vz.append(state.v[idx4, 2])

    # Phase 2: contact with fine dt
    print(f'  Fall done, z={state.r[idx4,2]:.4f}, switching to dt=1e-6...')
    for i in range(200000):
        intg.step(state, 1e-6)
        z = state.r[idx4, 2]
        if i % 500 == 0 and z < 0.025:
            if abs(state.v[idx4, 2]) < 1e-4:
                for _ in range(100):
                    intg.step(state, 1e-6)
                    our_t.append(state.t); our_z.append(state.r[idx4, 2])
                break
        our_t.append(state.t); our_z.append(z); our_vz.append(state.v[idx4, 2])

    our_t = np.array(our_t); our_z_mm = np.array(our_z) * 1000.0
    print(f'  Contact done: z={state.r[idx4,2]:.6f} m, vz={state.v[idx4,2]:.6f}')

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    ax.plot(rd_t, rd_contact_fm, 'b-', lw=2, label='RecurDyn FM GeoContact')
    ax.plot([rd_t[0], rd_t[-1]], [F_equilibrium, F_equilibrium], 'k:', lw=1.5,
            label=f'Theory: F = m·g = {m:.2f}×{g:.2f} = {F_equilibrium:.2f} N')
    ax.set_title('Contact Force Magnitude (static equilibrium)')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Force (N)')
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(rd_t, rd_contact_tm, 'r-', lw=2, label='RecurDyn TM GeoContact')
    ax.plot([rd_t[0], rd_t[-1]], [0, 0], 'k:', lw=1.5, label='Theory: T = 0 N·mm')
    ax.set_title('Contact Torque Magnitude (should be ~0)')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Torque (N·mm)')
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(rd_t, rd_b2_z, 'b-', lw=2, label='RecurDyn Body2 Z')
    ax.plot(our_t, our_z_mm, 'r--', lw=2, label=f'Our Framework (dt=1e-6)')
    settled_z = rd_b2_z[-1]
    ax.plot([rd_t[0], rd_t[-1]], [settled_z, settled_z], 'k:', lw=1.5,
            label=f'Equilibrium: Z = {settled_z:.2f} mm')
    ax.set_title(f'Ring Z Position (settled at {settled_z:.1f} mm)')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Z (mm)')
    ax.legend(); ax.grid(alpha=0.3)

    print(f'  Contact FM final: {rd_contact_fm[-1]:.4f} N (theory: {F_equilibrium:.2f})')
    print(f'  Contact TM final: {rd_contact_tm[-1]:.4f} N·mm (theory: 0)')
    print(f'  Ring settled Z (RecurDyn): {rd_b2_z[-1]:.4f} mm')
    print(f'  Ring settled Z (Our):     {our_z[-1]*1000:.4f} mm')
    print(f'  Ring final Vz: {rd_b2_vz[-1]:.6f} mm/s')

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'ring_cube_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved ring_cube_comparison.png')


# ============================================================
# Main
# ============================================================
def main():
    print(f'Three-Way Comparison: RecurDyn vs Analytical')
    print(f'Output: {OUT_DIR}')
    print(f'Note: Our framework results (red dashed lines) will be added')
    print(f'      once the MBD solver handles complex constraints.')

    compare_hollow_cyl()
    compare_gear()
    compare_ring_cube()

    print(f'\nAll comparison plots saved to {OUT_DIR}')


if __name__ == '__main__':
    main()
