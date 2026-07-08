import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.joints_revolute import revolute_joint_coordinate
from src.integrators.implicit_be import BackwardEulerIntegrator


def run():
    rmd_path = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'pendulum_real.rmd')

    converter = RMDToIR()
    model = converter.load_file(rmd_path)

    print(f"Case 08: RMD import (real format) + implicit Backward Euler")

    state = State(model)
    angle = 0.5
    c, s = np.cos(angle), np.sin(angle)
    state.R[1] = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
    state.r[1] = [0, 0, -np.cos(angle)]

    init_theta = revolute_joint_coordinate(model.joints[0], model, state)
    print(f"  Bodies: {model.num_bodies}, movable: {model.num_movable}")
    print(f"  Initial theta: {init_theta:.4f}")

    integrator = BackwardEulerIntegrator(model, tol=1e-6)
    T, dt = 0.5, 0.01
    history = integrator.integrate(state, T, dt)

    thetas = [r.get('joint_PinJoint_theta', 0) for r in history]
    print(f"  Steps: {len(history)}")
    print(f"  Theta range: [{min(thetas):.3f}, {max(thetas):.3f}]")
    print(f"  Final theta: {thetas[-1]:.4f}")
    assert len(history) == int(T / dt)
    print("  PASSED\n")
    return history


if __name__ == '__main__':
    run()
