import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import numpy as np
from src.ir import SystemModel, RigidBody, Frame, RevoluteJoint, FixedJoint
from src.ir import Drive, DriveMode, SinFunction, LinearFunction
from src.ir import ContactPair, ForceElement, ConstantFunction
from src.dynamics import State
from src.dynamics.joints_revolute import revolute_joint_coordinate
from src.integrators.explicit_projected import ExplicitProjectedIntegrator
from src.validation import RecurDynValidator, ValidationCase, compute_error_metrics


def run():
    print(f"Case 09: Framework validation against analytical solutions")
    print(f"{'='*60}")

    Rx = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    all_pass = True

    # ---- Test 1: Motion drive accuracy ----
    print(f"\n[Test 1] Motion drive: theta(t) = 2*PI*t")
    model1 = SystemModel()
    model1.add_body(RigidBody('g', 0, fixed=True))
    model1.add_body(RigidBody('p', 1, mass=1.0,
                              inertia_body=np.diag([0.1, 0.01, 0.1]),
                              initial_r=[0, 0, -1]))
    model1.add_frame(Frame('fg', 0, 0, local_orientation=Rx))
    model1.add_frame(Frame('fp', 1, 1, local_position=[0, 0, 1],
                           local_orientation=Rx))
    model1.add_joint(RevoluteJoint('pin', 0, 0, 1))
    model1.add_function('f', LinearFunction(0.0, 2*np.pi))
    model1.add_drive(Drive('d', 0, 0, mode=DriveMode.MOTION, function_id='f'))

    state1 = State(model1)
    intg1 = ExplicitProjectedIntegrator(model1)
    T1, dt1 = 0.2, 0.001
    h1 = intg1.integrate(state1, T1, dt1)

    thetas1 = [r.get('joint_pin_theta', 0) for r in h1]
    target1 = [2*np.pi*i*dt1 for i in range(len(thetas1))]
    err1 = max(abs(np.array(thetas1) - np.array(target1)))
    status1 = 'PASS' if err1 < 0.01 else 'FAIL'
    print(f"  Max tracking error: {err1:.6f} rad  [{status1}]")
    if status1 == 'FAIL': all_pass = False

    # ---- Test 2: External axial force ----
    print(f"\n[Test 2] Axial force: F = 10N, m = 1kg")
    model2 = SystemModel()
    model2.gravity = np.array([0, 0, 0])
    model2.add_body(RigidBody('g', 0, fixed=True))
    model2.add_body(RigidBody('b', 1, mass=1.0,
                              inertia_body=np.eye(3)*0.01,
                              initial_r=[0, 0, 0]))
    model2.add_frame(Frame('fi', 0, 1, local_position=[0, 0, 0]))
    model2.add_frame(Frame('fj', 1, 0, local_position=[0, 1, 0]))
    model2.add_force(ForceElement('axial', 0, 0, 1,
                                  ConstantFunction(10.0), 'axial'))

    state2 = State(model2)
    asm2 = __import__('src.dynamics.assembler', fromlist=['']).DynamicsAssembler(model2)
    _, _, Q2, _, _ = asm2.assemble(state2, 0.001)
    expected_F = 10.0  # 10N in Y direction
    actual_F = Q2[1]
    err2 = abs(actual_F - expected_F)
    status2 = 'PASS' if err2 < 1e-10 else 'FAIL'
    print(f"  Force on body: {actual_F:.2f} N (expected {expected_F:.2f} N)  [{status2}]")
    if status2 == 'FAIL': all_pass = False

    # ---- Test 3: Gravity loading ----
    print(f"\n[Test 3] Gravity: g = 9.81 m/s^2, m = 2kg")
    model3 = SystemModel()
    model3.gravity = np.array([0, 0, -9.81])
    model3.add_body(RigidBody('b', 0, mass=2.0, inertia_body=np.eye(3)*0.01))
    state3 = State(model3)
    asm3 = __import__('src.dynamics.assembler', fromlist=['']).DynamicsAssembler(model3)
    _, _, Q3, _, _ = asm3.assemble(state3, 0.001)
    expected_Fz = -2.0 * 9.81
    actual_Fz = Q3[2]
    err3 = abs(actual_Fz - expected_Fz)
    status3 = 'PASS' if err3 < 1e-10 else 'FAIL'
    print(f"  Gravity force: {actual_Fz:.2f} N (expected {expected_Fz:.2f} N)  [{status3}]")
    if status3 == 'FAIL': all_pass = False

    # ---- Test 4: Explicit integrator - free fall ----
    print(f"\n[Test 4] Free fall: z(t) = z0 - 0.5*g*t^2")
    model4 = SystemModel()
    model4.gravity = np.array([0, 0, -9.81])
    model4.add_body(RigidBody('b', 0, mass=1.0, inertia_body=np.eye(3)*0.01,
                              initial_r=[0, 0, 10]))
    state4 = State(model4)
    intg4 = ExplicitProjectedIntegrator(model4)
    h4 = intg4.integrate(state4, 1.0, 0.001)
    positions4 = [r.get('body_b_r', np.zeros(3))[2] for r in h4]
    t4 = np.array([i*0.001 for i in range(len(positions4))])
    expected_z = 10 - 0.5 * 9.81 * t4**2
    err4 = max(abs(np.array(positions4) - expected_z))
    status4 = 'PASS' if err4 < 0.05 else 'FAIL'
    print(f"  Max position error: {err4:.4f} m  [{status4}]")
    if status4 == 'FAIL': all_pass = False

    # ---- Test 5: RMD import - read real format ----
    print(f"\n[Test 5] RMD real format import")
    rmd_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
    rmd_path = os.path.join(rmd_dir, 'pendulum_real.rmd')
    from src.importers import RMDToIR
    converter = RMDToIR()
    model5 = converter.load_file(rmd_path)
    status5 = 'PASS' if model5.num_bodies == 2 and model5.num_movable == 1 else 'FAIL'
    print(f"  Bodies={model5.num_bodies} movable={model5.num_movable}  [{status5}]")
    if status5 == 'FAIL': all_pass = False

    print(f"\n{'='*60}")
    print(f"  Overall: {'ALL PASSED' if all_pass else 'SOME FAILED'}")
    print(f"{'='*60}")
    assert all_pass, "Validation tests failed"
    print("  PASSED\n")
    return all_pass


if __name__ == '__main__':
    run()
