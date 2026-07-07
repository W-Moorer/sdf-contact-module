#!/usr/bin/env python3
"""Level 0: Analytic plane SDF validation — ring-on-cube torsion.

Runs analytic quadrature + all OBJ mesh variants and generates
convergence plot.
"""

import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from src.scenes.ring_on_cube_torsion import run_all_comparisons
from src.visualization import save_convergence_plot

if __name__ == '__main__':
    results = run_all_comparisons('outputs')

    print("\nGenerating convergence plot...")
    save_convergence_plot(results, 'outputs/ring_cube_level0')

    print("\nAll outputs saved to outputs/")
