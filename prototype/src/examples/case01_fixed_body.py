import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import SystemModel, RigidBody, Frame, FixedJoint
from src.dynamics import State
from src.integrators import ExplicitProjectedIntegrator


def build_model():
    model = SystemModel()

    body = RigidBody(name='block', body_id=0, mass=1.0,
                     inertia_body=np.eye(3) * 0.1,
                     initial_r=[0, 0, 0],
                     fixed=True)
    model.add_body(body)

    frame_w = Frame(name='world', frame_id=0, body_id=0,
                    local_position=[0, 0, 0])
    frame_b = Frame(name='block_frame', frame_id=1, body_id=0,
                    local_position=[0, 0, 0])
    model.add_frame(frame_w)
    model.add_frame(frame_b)

    j = FixedJoint(name='fix_block', joint_id=0,
                   frame_i_id=0, frame_j_id=1)
    model.add_joint(j)

    return model


def run():
    model = build_model()
    state = State(model)
    integrator = ExplicitProjectedIntegrator(model)

    T = 1.0
    dt = 0.01
    integrator.integrate(state, T, dt)

    print(f"Case 01: Fixed body (static check)")
    print(f"  Final position: {state.r[0]}")
    print(f"  Position error: {np.linalg.norm(state.r[0]):.2e}")
    assert np.linalg.norm(state.r[0]) < 1e-6, "Fixed body should not move"
    print("  PASSED\n")
    return integrator.history


if __name__ == '__main__':
    run()
