#!/usr/bin/env python3
"""Level 0: Analytic plane SDF validation — ring-on-cube torsion."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from src.scenes.ring_on_cube_torsion import run_ring_on_cube_torsion

if __name__ == '__main__':
    diag = run_ring_on_cube_torsion(output_dir='outputs/ring_cube_level0')
    print("\nAll figures saved to outputs/ring_cube_level0/figures/")
