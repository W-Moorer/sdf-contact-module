"""Granular profile of assemble + projection internals."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR
from src.dynamics import State
from src.dynamics.assembler import DynamicsAssembler
from src.dynamics.constraints import (_all_residuals, num_joint_constraints,
                                      num_drive_constraints, eval_constraint_jacobian)
from src.dynamics.drives import motion_drive_jacobian
from src.dynamics.mass_matrix import build_mass_matrix
from src.integrators.projection import _build_jacobian_full, _pseudo_inverse_solve
from src.solvers.kkt import solve_kkt

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
                    'Ring - Cube - Centric Collision_01',
                    'Ring - Cube - Centric Collision.rmd')

converter = RMDToIR()
model = converter.load_file(RMD)
state = State(model)
asm = DynamicsAssembler(model)

nc = num_joint_constraints(model) + num_drive_constraints(model)
nv = 6 * model.num_movable
print(f'Bodies: {model.num_movable}, nc={nc}, nv={nv}')

import time as _time
times = {
    'build_mass_matrix':[], 'eval_J':[], 'build_J':[], 'motion_drive_J':[],
    'compute_bc':[], 'J_to_sparse':[],
    'kkt':[], 'integrate':[],
    'proj_pos_total':[], 'proj_pos_niter':[],
    'proj_vel_build_J':[], 'proj_vel_mass':[], 'proj_vel_svd':[], 'proj_vel_total':[],
    'proj_pos_residual':[], 'proj_pos_build_J':[], 'proj_pos_mass':[], 'proj_pos_svd':[]
}

for step in range(100):
    # ---- ASSEMBLE breakdown ----
    t0 = _time.perf_counter()
    M = build_mass_matrix(model, state)
    times['build_mass_matrix'].append(_time.perf_counter() - t0)

    t0 = _time.perf_counter()
    J_joint = eval_constraint_jacobian(model, state)
    times['eval_J'].append(_time.perf_counter() - t0)

    t0 = _time.perf_counter()
    nd = num_drive_constraints(model)
    if nd > 0:
        J_drive = np.zeros((nd, nv))
        row = 0
        for drive in model.drives.values():
            Jd = motion_drive_jacobian(drive, model, state)
            J_drive[row:row+1] = Jd
            row += 1
        J = np.vstack([J_joint, J_drive]) if len(J_joint) > 0 else J_drive
    else:
        J = J_joint
    times['motion_drive_J'].append(_time.perf_counter() - t0)

    t0 = _time.perf_counter()
    from src.dynamics.constraints import compute_bc_from_jacobian
    b_c = compute_bc_from_jacobian(J, model, state, 1e-6)
    times['compute_bc'].append(_time.perf_counter() - t0)

    from scipy.sparse import csr_matrix
    t0 = _time.perf_counter()
    Jsp = csr_matrix(J)
    Msp = csr_matrix(M)
    times['J_to_sparse'].append(_time.perf_counter() - t0)

    # ---- KKT ----
    t0 = _time.perf_counter()
    C = np.zeros(nv)  # gyroscopic is zero for this model
    Q_ext = asm._external_forces(state)
    acc, lam = solve_kkt(Msp, C, Q_ext, Jsp, b_c)
    times['kkt'].append(_time.perf_counter() - t0)

    # ---- integrate ----
    t0 = _time.perf_counter()
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
            c = np.cos(dtheta_norm); s = np.sin(dtheta_norm)
            dR = (c*np.eye(3) + s*np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]) + (1-c)*np.outer(axis,axis))
            state.R[idx] = dR @ state.R[idx]
        tmp_idx += 1
    times['integrate'].append(_time.perf_counter() - t0)

    # ---- POSITION PROJECTION ----
    t0 = _time.perf_counter()
    for iteration in range(20):
        Phi = _all_residuals(model, state)
        err = np.max(np.abs(Phi))
        if err < 1e-10: break
        t1 = _time.perf_counter()
        Jp = _build_jacobian_full(model, state)
        times.setdefault('proj_pos_build_J',[]).append(_time.perf_counter()-t1)
        t1 = _time.perf_counter()
        Mp = build_mass_matrix(model, state).toarray()
        times.setdefault('proj_pos_mass',[]).append(_time.perf_counter()-t1)
        t1 = _time.perf_counter()
        M_inv = np.linalg.inv(Mp)
        W = Jp @ M_inv @ Jp.T
        lam_p = _pseudo_inverse_solve(W, -Phi)
        delta_V = M_inv @ Jp.T @ lam_p
        times.setdefault('proj_pos_svd',[]).append(_time.perf_counter()-t1)
        tmp_idx = 0
        for body in model.movable_bodies:
            idx = state.body_idx(body.id)
            state.r[idx] += delta_V[6*tmp_idx:6*tmp_idx+3]
            dtheta = delta_V[6*tmp_idx+3:6*tmp_idx+6]
            dtheta_norm = np.linalg.norm(dtheta)
            if dtheta_norm > 1e-30:
                axis = dtheta / dtheta_norm
                c,s = np.cos(dtheta_norm), np.sin(dtheta_norm)
                dR = c*np.eye(3) + s*np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]]) + (1-c)*np.outer(axis,axis)
                state.R[idx] = dR @ state.R[idx]
            tmp_idx += 1
    times['proj_pos_total'].append(_time.perf_counter() - t0)
    times['proj_pos_niter'].append(iteration)

    # ---- VELOCITY PROJECTION ----
    t0 = _time.perf_counter()
    t1 = _time.perf_counter()
    Jv = _build_jacobian_full(model, state)
    times['proj_vel_build_J'].append(_time.perf_counter()-t1)
    t1 = _time.perf_counter()
    Mv = build_mass_matrix(model, state).toarray()
    times['proj_vel_mass'].append(_time.perf_counter()-t1)
    t1 = _time.perf_counter()
    V = state.pack_V()
    M_inv_v = np.linalg.inv(Mv)
    lam_v = _pseudo_inverse_solve(Jv @ M_inv_v @ Jv.T, -(Jv @ V))
    state.unpack_V(V + M_inv_v @ Jv.T @ lam_v)
    times['proj_vel_svd'].append(_time.perf_counter()-t1)
    times['proj_vel_total'].append(_time.perf_counter() - t0)

    state.t += 1e-6

print(f"\n{'section':>24s} {'avg(us)':>10s} {'total(ms)':>12s} {'calls':>6s}")
print("-"*55)
for k, lbl in [
    ('build_mass_matrix','build_mass_matrix'), ('eval_J','eval_constraint_J'),
    ('motion_drive_J','motion_drive_J'),
    ('compute_bc','compute_bc'), ('J_to_sparse','J→sparse'),
    ('kkt','kkt'), ('integrate','integrate'),
    ('proj_pos_total','pos_proj (total)'), ('proj_pos_niter','  pos_proj iter'),
    ('proj_pos_build_J','  pos build_J'),
    ('proj_pos_mass','  pos mass'),
    ('proj_pos_svd','  pos SVD'),
    ('proj_vel_total','vel_proj (total)'),
    ('proj_vel_build_J','  vel build_J'),
    ('proj_vel_mass','  vel mass'),
    ('proj_vel_svd','  vel SVD'),
]:
    arr = np.array(times[k])
    avg_us = np.mean(arr) * 1e6
    total_ms = np.sum(arr) * 1e3
    n = len(arr)
    print(f"{lbl:>24s} {avg_us:>10.1f} {total_ms:>12.2f} {n:>6d}")

tot = sum(np.sum(np.array(v)) for k,v in times.items() if k not in ['proj_pos_niter'])
print(f"\nTotal avg/step: {tot/100*1e6:.0f} us")
print(f"Theoretical rate: {1/(tot/100):.0f} steps/s")
