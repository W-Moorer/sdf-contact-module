import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import (SystemModel, RigidBody, Frame, RevoluteJoint,
                    Drive, DriveMode, SinFunction)
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

    model.add_function('drv', SinFunction(amplitude=np.pi/4, omega=2.0))
    model.add_drive(Drive('swing', 0, target_joint_id=0,
                          mode=DriveMode.MOTION, function_id='drv'))
    return model


def run():
    model = build_model()
    state = State(model)
    integrator = ExplicitProjectedIntegrator(model)

    history = integrator.integrate(state, 3.0, 0.001)

    thetas = [r.get('joint_pin_theta', 0) for r in history]
    print(f"Case 03: Driven revolute (sinusoidal drive)")
    print(f"  Steps: {len(history)}")
    print(f"  Range: [{min(thetas):.3f}, {max(thetas):.3f}]")
    assert max(abs(a) for a in thetas) > 0.1
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
