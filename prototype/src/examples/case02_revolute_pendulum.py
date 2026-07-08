import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import SystemModel, RigidBody, Frame, RevoluteJoint
from src.dynamics import State
from src.integrators import ExplicitProjectedIntegrator


def build_model():
    model = SystemModel()
    L = 1.0
    model.add_body(RigidBody('ground', 0, fixed=True))
    model.add_body(RigidBody('pendulum', 1, mass=1.0,
                             inertia_body=np.diag([0.01, 0.1, 0.01]),
                             initial_r=[0, 0, -L]))
    Rx = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    model.add_frame(Frame('fg', 0, 0, local_orientation=Rx))
    model.add_frame(Frame('fp', 1, 1, local_position=[0, 0, L],
                          local_orientation=Rx))
    model.add_joint(RevoluteJoint('pin', 0, 0, 1))
    return model


def run():
    model = build_model()
    state = State(model)

    angle = 0.5; c, s = np.cos(angle), np.sin(angle)
    state.R[1] = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
    state.r[1] = [-s, 0, -c]

    integrator = ExplicitProjectedIntegrator(model)
    history = integrator.integrate(state, 3.0, 0.001)

    thetas = [r.get('joint_pin_theta', 0) for r in history]
    print(f"Case 02: Revolute pendulum (swinging)")
    print(f"  Steps: {len(history)}, Init: {thetas[0]:.3f}")
    print(f"  Range: [{min(thetas):.3f}, {max(thetas):.3f}]")
    assert max(abs(a) for a in thetas) > 0.4
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
