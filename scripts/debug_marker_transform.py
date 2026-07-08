"""Debug the marker-based transform at contact height."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers import RMDToIR; from src.dynamics import State
from src.dynamics.contacts_sdf import SDFContactEngine
from src.geometry.quadrature import QuadratureMesh
from src.sdf_grid import TrilinearSDFGrid

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')
converter = RMDToIR(); model = converter.load_file(RMD)
for cp in model.contacts.values():
    cp.body_a_id, cp.body_b_id = cp.body_b_id, cp.body_a_id
    cp.action_marker_id, cp.base_marker_id = cp.base_marker_id, cp.action_marker_id

cp = list(model.contacts.values())[0]
print(f'action_marker_id={cp.action_marker_id} base_marker_id={cp.base_marker_id}')
af = model.frames.get(cp.action_marker_id)
bf = model.frames.get(cp.base_marker_id)
print(f'Action: "{af.name}" body={af.body_id} pos={af.local_position}')
print(f'Base:   "{bf.name}" body={bf.body_id} pos={bf.local_position}')

state = State(model); idx4 = state.body_idx(4)
for ring_z in [0.250, 0.100, 0.065, 0.060, 0.000]:
    state.r[idx4, 2] = ring_z
    R, t = af.relative_to(bf, state)
    z_base = t[2] + 0.009  # lowest quad point z in base frame
    gap = z_base - 0.13     # cube top at z=0.13 in base frame
    print(f'ring_z={ring_z*1000:.0f}mm -> z_base={z_base*1000:.1f}mm gap={gap*1000:.1f}mm')
