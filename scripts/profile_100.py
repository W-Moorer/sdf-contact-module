"""Profile 100 steps to find bottleneck."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.assembler import DynamicsAssembler
from src.solvers.kkt import solve_kkt
from src.integrators.projection import position_projection, velocity_projection

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
                    'Ring - Cube - Centric Collision_01',
                    'Ring - Cube - Centric Collision.rmd')

converter = RMDToIR()
model = converter.load_file(RMD)
state = State(model)
asm = DynamicsAssembler(model)

times = {'assemble':[], 'kkt':[], 'integrate':[], 'proj_pos':[], 'proj_vel':[]}
for _ in range(100):
    t0 = time.perf_counter()
    M, C, Q, J, b_c = asm.assemble(state, 1e-6)
    times['assemble'].append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    acc, lam = solve_kkt(M, C, Q, J, b_c)
    times['kkt'].append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    tmp_idx = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        state.v[idx] += 1e-6 * acc[6*tmp_idx:6*tmp_idx+3]
        state.omega[idx] += 1e-6 * acc[6*tmp_idx+3:6*tmp_idx+6]
        tmp_idx += 1
    tmp_idx = 0
    for body in model.movable_bodies:
        idx = state.body_idx(body.id)
        state.r[idx] += 1e-6 * state.v[idx]
        dtheta = 1e-6 * state.omega[idx]
        dtheta_norm = np.linalg.norm(dtheta)
        if dtheta_norm > 1e-30:
            axis = dtheta / dtheta_norm
            c = np.cos(dtheta_norm)
            s = np.sin(dtheta_norm)
            dR = (c * np.eye(3) + s * np.array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],[-axis[1], axis[0], 0]]) + (1 - c) * np.outer(axis, axis))
            state.R[idx] = dR @ state.R[idx]
        tmp_idx += 1
    times['integrate'].append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    position_projection(model, state)
    times['proj_pos'].append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    velocity_projection(model, state)
    times['proj_vel'].append(time.perf_counter() - t0)

    state.t += 1e-6

print(f"{'step':>12s} {'avg(us)':>10s} {'total(ms)':>10s} {'%':>6s}")
print("-"*40)
for k in times:
    arr = np.array(times[k])
    avg_us = np.mean(arr) * 1e6
    total_ms = np.sum(arr) * 1e3
    pct = np.sum(arr) / sum(np.sum(np.array(v)) for v in times.values()) * 100
    print(f"{k:>12s} {avg_us:>10.1f} {total_ms:>10.2f} {pct:>5.1f}%")
print(f"\nTotal avg/step: {sum(np.mean(np.array(v)) for v in times.values())*1e6:.1f} us")
print(f"Theoretical rate: {1e6/sum(np.mean(np.array(v)) for v in times.values())*1e6:.0f} steps/s")
