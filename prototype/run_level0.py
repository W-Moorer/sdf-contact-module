#!/usr/bin/env python3
"""Level 0–1 validation: analytic plane SDF, OBJ mesh convergence,
and trilinear SDF grid convergence."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from src.scenes.ring_on_cube_torsion import run_all

if __name__ == '__main__':
    run_all('outputs')
    print("\nAll outputs saved to outputs/")
