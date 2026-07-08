"""Verify GSurface coordinate transformation from RM_MARKER to body COM."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1, errors='replace')
import numpy as np
from src.importers.rmd_lexer import RMDLexer

MODELS = os.path.join(os.path.dirname(__file__), '..', 'models')
RMD = os.path.join(MODELS, '圆环-立方体-对心碰撞',
    'Ring - Cube - Centric Collision_01', 'Ring - Cube - Centric Collision.rmd')

lexer = RMDLexer()
entities = lexer.load(RMD)

# --- 1. Find all markers and their parts ---
print("=== MARKERS ===")
for ent in entities:
    if ent.kind == 'MARKER':
        attrs = ent.attrs
        name = attrs.get('name', f"MARKER_{ent.id}")
        part_id = attrs.get('part_id', attrs.get('pid', '?'))
        pos = (attrs.get('xp', attrs.get('position', None)) or 
               attrs.get('x1', attrs.get('orig_x', [0,0,0])))
        ori_euler = attrs.get('euler', None)
        print(f"  ID={ent.id} '{name}' part={part_id} pos={pos} ori={ori_euler}")

# --- 2. Find GGEOM entities and their RM_MARKER ---
print("\n=== GGEOM (Surfaces) ===")
ggeom_marker_map = {}
for ent in entities:
    if ent.kind == 'GGEOM':
        attrs = ent.attrs
        name = attrs.get('name', f"GGEOM_{ent.id}")
        rm_marker = attrs.get('rm_marker', attrs.get('rm_id', attrs.get('marker', None)))
        print(f"  ID={ent.id} '{name}' RM_MARKER={rm_marker}")
        
        # Parse NODES to find z-range
        nodes_z = []
        mode = None
        for line in ent.raw_lines:
            s = line.strip()
            if not s or s.startswith('!'): continue
            if ', NODES' in s.upper(): mode = 'nodes'; continue
            if ', PATCHES' in s.upper(): mode = 'patches'; continue
            if mode == 'nodes':
                parts = [p.strip() for p in s.lstrip(',').split(',') if p.strip()]
                if len(parts) == 3:
                    try:
                        x, y, z = map(float, parts)
                        nodes_z.append(z)
                    except: pass
        
        if nodes_z:
            print(f"    NODES z-range: [{min(nodes_z):.4f}, {max(nodes_z):.4f}]")
        ggeom_marker_map[ent.id] = rm_marker

# --- 3. Build the body → marker mapping ---
print("\n=== PARTS ===")
for ent in entities:
    if ent.kind == 'PART':
        print(f"  ID={ent.id} name={ent.attrs.get('name','?')}")

# --- 4. Check: which marker is the RM_MARKER for the ring's surface? ---
# From RMD: Body2 (ring, id=4) has the contact surface
# The GGEOM with "Ring" or similar name
print("\n=== CONTACT SURFACE ANALYSIS ===")
from src.importers import RMDToIR
converter = RMDToIR()
model = converter.load_file(RMD)

# Check the contact pair body IDs
for cid, cp in model.contacts.items():
    ba = model.bodies[cp.body_a_id]
    bb = model.bodies[cp.body_b_id]
    print(f"Contact {cid}: body_a={ba.name}(id={ba.id}) body_b={bb.name}(id={bb.id})")

# Check frames attached to the ring body (body 4)
print("\n=== FRAMES (Markers) attached to ring body ===")
for fid, frame in model.frames.items():
    if frame.body_id == 4:  # ring body
        print(f"  Frame {fid}: name={frame.name} pos={frame.local_position} ori={frame.local_orientation}")

# Check frames attached to the cube body (body 3)
print("\n=== FRAMES (Markers) attached to cube body ===")
for fid, frame in model.frames.items():
    if frame.body_id == 3:  # cube body
        print(f"  Frame {fid}: name={frame.name} pos={frame.local_position} ori={frame.local_orientation}")
