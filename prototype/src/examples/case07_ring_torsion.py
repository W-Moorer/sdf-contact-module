import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import SystemModel, RigidBody, ContactPair
from src.dynamics import State
from src.integrators import ExplicitProjectedIntegrator
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid
from src.analytic_benchmarks import analytic_torsion_torque


def build_annular_quadrature(R_i, R_o, n_radial, n_azimuthal):
    dr = (R_o - R_i) / n_radial
    dtheta = 2.0 * np.pi / n_azimuthal
    X_q, w_q = [], []
    for i in range(n_radial):
        r_inner = R_i + i * dr
        r_outer = R_i + (i + 1) * dr
        r_center = 0.5 * (r_inner + r_outer)
        area = 0.5 * (r_outer**2 - r_inner**2) * dtheta
        for j in range(n_azimuthal):
            theta = (j + 0.5) * dtheta
            X_q.append([r_center * np.cos(theta), r_center * np.sin(theta), 0.0])
            w_q.append(area)
    return QuadratureMesh(np.array(X_q), np.array(w_q))


def run():
    model = SystemModel()
    model.gravity = np.array([0, 0, -9.81])

    cube = RigidBody(name='cube', body_id=0, fixed=True)
    model.add_body(cube)

    R_i, R_o = 0.4, 0.6
    mass = 10.0
    cylinder = RigidBody(name='cylinder', body_id=1, mass=mass,
                         inertia_body=np.diag([0.1, 0.1, 0.1]),
                         initial_r=[0, 0, 0.05],
                         initial_omega=[0, 0, 1.0])
    model.add_body(cylinder)

    cp = ContactPair(name='ring_contact', contact_id=0,
                     body_a_id=1, body_b_id=0,
                     normal_stiffness=1e6, friction_coefficient=0.5,
                     activation_distance=0.0,
                     quadrature_settings={'regularizer': 1e-4})
    model.add_contact(cp)

    quad_mesh = build_annular_quadrature(R_i, R_o, 16, 48)

    models_dir = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'models')
    sdf_grid = TrilinearSDFGrid(os.path.join(models_dir, 'cube_res128.sdf'))

    engine = SDFContactEngine()
    engine.register_pair(0, quad_mesh, sdf_grid)

    integrator = ExplicitProjectedIntegrator(model, contact_engine=engine)
    state = State(model)

    T, dt = 0.2, 0.0005
    history = integrator.integrate(state, T, dt)

    positions = [r.get('body_cylinder_r', np.zeros(3))[2] for r in history]
    print(f"Case 07: Ring-on-cube torsion (trilinear SDF)")
    print(f"  Steps: {len(history)}")
    print(f"  Initial z: {positions[0]:.4f}")
    print(f"  Final z:   {positions[-1]:.4f}")
    print(f"  Min z:     {min(positions):.4f}")

    steady_z = np.mean(positions[-50:])
    print(f"  Settled z: {steady_z:.4f} (penetration: {-steady_z:.4f})")

    assert min(positions) > -0.02, f"Penetrated too deep: {min(positions)}"
    assert len(history) > 100
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
