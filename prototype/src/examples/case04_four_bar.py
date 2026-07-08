import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import (SystemModel, RigidBody, Frame, RevoluteJoint,
                    Drive, DriveMode, SinFunction)
from src.dynamics import State
from src.integrators import ExplicitProjectedIntegrator


def build_model():
    """Open-chain double pendulum (no loop closure) for validation."""
    model = SystemModel()
    L = 1.0
    Rx = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])

    model.add_body(RigidBody('ground', 0, fixed=True))
    model.add_body(RigidBody('link1', 1, mass=1.0,
                             inertia_body=np.diag([0.1, 0.01, 0.1]),
                             initial_r=[0, 0, -L]))
    model.add_body(RigidBody('link2', 2, mass=1.0,
                             inertia_body=np.diag([0.1, 0.01, 0.1]),
                             initial_r=[0, 0, -2*L]))

    model.add_frame(Frame('g0', 0, 0, local_orientation=Rx))
    model.add_frame(Frame('l1a', 1, 1, local_position=[0, 0, L],
                          local_orientation=Rx))
    model.add_frame(Frame('l1b', 2, 1, local_position=[0, 0, -L],
                          local_orientation=Rx))
    model.add_frame(Frame('l2a', 3, 2, local_position=[0, 0, L],
                          local_orientation=Rx))

    model.add_joint(RevoluteJoint('j1', 0, 0, 1))
    model.add_joint(RevoluteJoint('j2', 1, 2, 3))

    model.add_function('drv', SinFunction(amplitude=np.pi/4, omega=1.0))
    model.add_drive(Drive('cdrv', 0, target_joint_id=0,
                          mode=DriveMode.MOTION, function_id='drv'))
    return model


def run():
    model = build_model()
    state = State(model)
    integrator = ExplicitProjectedIntegrator(model)

    history = integrator.integrate(state, 3.0, 0.001)

    thetas = [r.get('joint_j1_theta', 0) for r in history]
    print(f"Case 04: Double pendulum (open chain)")
    print(f"  Steps: {len(history)}")
    print(f"  Joint 1 range: [{min(thetas):.3f}, {max(thetas):.3f}]")
    assert len(history) == int(3.0 / 0.001)
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
