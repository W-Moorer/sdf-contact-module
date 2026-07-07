#!/usr/bin/env python3
"""Verify hollow_cylinder SDF accuracy at different resolutions.

For each surface triangle centroid, query the SDF and compare:
  - gap vs 0 (should be 0 on the surface)
  - gradient direction vs mesh face normal
"""

import numpy as np
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from src.sdf_grid import TrilinearSDFGrid
from src.mesh_quadrature import load_obj_triangles, triangle_centroid_quadrature

models_dir = os.path.join(os.path.dirname(__file__), '..', 'models')

def check_sdf_accuracy(obj_path, sdf_path, label):
    verts, tris = load_obj_triangles(obj_path)
    X_q, w_q, tri_normals, _ = triangle_centroid_quadrature(verts, tris)

    grid = TrilinearSDFGrid(sdf_path)

    gaps = np.zeros(len(X_q))
    grad_dot_n = np.zeros(len(X_q))
    grad_norms = np.zeros(len(X_q))

    for i, p in enumerate(X_q):
        g, grad, gn, n, valid = grid.query(p)
        gaps[i] = g
        grad_norms[i] = gn
        if gn > 1e-30:
            unit_grad = grad / gn
        else:
            unit_grad = np.zeros(3)
        grad_dot_n[i] = np.dot(unit_grad, tri_normals[i])

    gap_rms = np.sqrt(np.mean(gaps**2))
    gap_max = np.max(np.abs(gaps))
    dot_mean = np.mean(grad_dot_n)
    dot_std = np.std(grad_dot_n)
    gn_mean = np.mean(grad_norms)
    gn_std = np.std(grad_norms)

    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"{'='*60}")
    print(f"  SDF: {os.path.basename(sdf_path)}")
    print(f"  Points: {len(X_q)}")
    print(f"  Gap RMS:         {gap_rms:.6e}")
    print(f"  Gap |max|:        {gap_max:.6e}")
    print(f"  grad·n_face mean: {dot_mean:.6f}  (should be ≈1)")
    print(f"  grad·n_face std:  {dot_std:.6f}")
    print(f"  |grad| mean:      {gn_mean:.6f}  (should be ≈1)")
    print(f"  |grad| std:       {gn_std:.6f}")

    # Categorize by surface type
    tol = 0.01
    bottom_mask = tri_normals[:, 2] < -(1 - tol)
    top_mask = tri_normals[:, 2] > (1 - tol)
    outer_mask = (~bottom_mask) & (~top_mask)

    print(f"\n  --- By surface type ---")
    for name, mask in [("bottom face", bottom_mask), ("top face", top_mask),
                       ("cylindrical wall", outer_mask)]:
        if mask.sum() == 0:
            continue
        g_rms = np.sqrt(np.mean(gaps[mask]**2))
        g_max = np.max(np.abs(gaps[mask]))
        d = np.mean(grad_dot_n[mask])
        gf = np.mean(grad_norms[mask])
        print(f"  {name:20s}:  gap_rms={g_rms:.2e}  gap_max={g_max:.2e}  "
              f"grad·n={d:.4f}  |grad|={gf:.4f}")

    return {
        'label': label,
        'sdf_file': os.path.basename(sdf_path),
        'num_points': len(X_q),
        'gap_rms': float(gap_rms),
        'gap_max': float(gap_max),
        'grad_dot_n_mean': float(dot_mean),
        'grad_dot_n_std': float(dot_std),
        'grad_norm_mean': float(gn_mean),
        'grad_norm_std': float(gn_std),
    }

if __name__ == '__main__':
    results = []
    for res in [64, 128, 256]:
        sdf_path = os.path.join(models_dir, f"cylinder_res{res}.sdf")
        if os.path.exists(sdf_path):
            d = check_sdf_accuracy(
                os.path.join(models_dir, 'hollow_cylinder_fine.obj'),
                sdf_path, f"Resolution {res}")
            results.append(d)

    print(f"\n{'='*60}")
    print(f"Summary")
    print(f"{'='*60}")
    print(f"{'Res':>6s}  {'Gap RMS':>10s}  {'Gap |max|':>10s}  "
          f"{'grad·n':>8s}  {'|grad|':>8s}")
    for r in results:
        print(f"{r['label'][-3:]:>6s}  {r['gap_rms']:>10.2e}  "
              f"{r['gap_max']:>10.2e}  {r['grad_dot_n_mean']:>8.4f}  "
              f"{r['grad_norm_mean']:>8.4f}")

    # Also do the full Level 2: cylinder OBJ integration against cube SDF
    print(f"\n{'='*60}")
    print(f"Level 2: cylinder OBJ vs cube SDF (contact integration)")
    print(f"{'='*60}")
    from src.scenes.ring_on_cube_torsion import run_level1b_obj, run_level0
    print(f"\nReference (analytic plane):")
    ref = run_level0('outputs/check_cylinder')
    print(f"\nTrilinear cube SDF + cylinder OBJ:")
    for sdf_res in [64, 128, 256]:
        sdf_path = os.path.join(models_dir, f"cube_res{sdf_res}.sdf")
        if os.path.exists(sdf_path):
            d = run_level1b_obj('outputs/check_cylinder',
                os.path.join(models_dir, 'hollow_cylinder_fine.obj'),
                sdf_path, 'fine', sdf_res)

    print("\nDone.")
