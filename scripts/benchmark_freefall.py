"""Benchmark: 1M steps dt=1e-6, NO contact, just free fall."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.integrators.explicit_projected import ExplicitProjectedIntegrator

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
                    'Ring - Cube - Centric Collision_01',
                    'Ring - Cube - Centric Collision.rmd')

converter = RMDToIR()
model = converter.load_file(RMD)
state = State(model)
intg = ExplicitProjectedIntegrator(model)

N = 1_000_000
t0 = time.perf_counter()
for i in range(N):
    intg.step(state, 1e-6)
    if i % 200_000 == 0 and i > 0:
        elapsed = time.perf_counter() - t0
        rate = i / elapsed
        print(f'  {i//1000:4d}k t={state.t:.4f}s rate={rate:.0f} steps/s')
elapsed = time.perf_counter() - t0
print(f'\n{N} steps in {elapsed:.1f}s => {N/elapsed:.0f} steps/s, {elapsed/N*1e6:.1f} us/step')
print(f'Final: t={state.t:.6f}, z={state.r[state.body_idx(4),2]:.6f}')
