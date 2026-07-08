import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import SystemModel, RigidBody, ContactPair
from src.dynamics import State
from src.integrators import ExplicitProjectedIntegrator
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import AnalyticPlaneSDF


def run():
    model = SystemModel()
    model.gravity = np.array([0, 0, -9.81])

    ground = RigidBody(name='ground', body_id=0, fixed=True)
    model.add_body(ground)

    block = RigidBody(name='block', body_id=1, mass=1.0,
                      inertia_body=np.eye(3) * 0.01,
                      initial_r=[0, 0, 0.05])
    model.add_body(block)

    cp = ContactPair(name='contact', contact_id=0,
                     body_a_id=1, body_b_id=0,
                     normal_stiffness=1e6, friction_coefficient=0.5,
                     activation_distance=0.0,
                     quadrature_settings={'regularizer': 1e-4})
    model.add_contact(cp)

    side = 0.3
    n_grid = 3
    xs = np.linspace(-side/2, side/2, n_grid)
    X_q = np.array([[x, y, 0.0] for x in xs for y in xs])
    w_q = np.full(len(X_q), (side / n_grid) ** 2)
    quad_mesh = QuadratureMesh(X_q, w_q)

    sdf_ground = AnalyticPlaneSDF(z_top=0.0)

    engine = SDFContactEngine()
    engine.register_pair(0, quad_mesh, sdf_ground)

    integrator = ExplicitProjectedIntegrator(model, contact_engine=engine)
    state = State(model)

    T, dt = 0.15, 0.0005
    history = integrator.integrate(state, T, dt)

    positions = [r.get('body_block_r', np.zeros(3))[2] for r in history]
    print(f"Case 06: Block falling onto plane (SDF contact in MBD loop)")
    print(f"  Steps: {len(history)}")
    print(f"  Initial z: {positions[0]:.4f}")
    print(f"  Final z:   {positions[-1]:.4f}")
    print(f"  Min z:     {min(positions):.4f}")
    print(f"  Settled: {abs(positions[-1]) < 0.02}")

    assert min(positions) > -0.01, f"Block penetrated too deep: {min(positions)}"
    assert len(history) > 100
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
