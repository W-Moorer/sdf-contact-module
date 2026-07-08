import os
import numpy as np


def parse_patch_line(line):
    parts = line.strip().split(',')
    nums = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        try:
            nums.append(float(p))
        except ValueError:
            pass
    return nums


def extract_mesh_from_rmd(entities, ggeom_id):
    ggeoms = [e for e in entities if e.kind == 'GGEOM' and e.id == ggeom_id]
    if not ggeoms:
        return None
    ent = ggeoms[0]

    npatch = ent.get_int('NO_PATCH')
    nnodes = ent.get_int('NO_NODE')
    ncell = ent.get_floats('NCELL')
    cell_size = ent.get_floats('CELL_SIZE')

    in_patches = False
    vertices = {}
    triangles = []
    vertex_counter = 1

    for line in ent.raw_lines:
        s = line.strip()
        if 'PATCHES' in s.upper():
            in_patches = True
            continue
        if not in_patches:
            continue

        nums = parse_patch_line(s)
        if len(nums) < 4:
            continue

        patch_type = int(nums[0])
        if patch_type != 3:
            continue

        n1, n2, n3 = int(nums[1]), int(nums[2]), int(nums[3])

        for nid in (n1, n2, n3):
            if nid not in vertices:
                vertices[nid] = vertex_counter
                vertex_counter += 1

        triangles.append((vertices[n1], vertices[n2], vertices[n3]))

    if len(triangles) == 0:
        return None

    return {
        'num_patches': npatch,
        'num_nodes': nnodes,
        'num_triangles': len(triangles),
        'num_unique_vertices': len(vertices),
        'vertices': vertices,
        'triangles': triangles,
        'ncell': ncell,
        'cell_size': cell_size,
    }


def assign_grid_positions(mesh_data, length_scale=1.0):
    ncell = mesh_data['ncell']
    cell_size = mesh_data['cell_size']
    nnodes = mesh_data['num_nodes']

    if len(ncell) < 3 or len(cell_size) < 3:
        return np.zeros((len(mesh_data['vertices']), 3))

    nx, ny, nz = int(ncell[0]), int(ncell[1]), int(ncell[2])
    dx = cell_size[0] / nx if nx > 0 else 1.0
    dy = cell_size[1] / ny if ny > 0 else 1.0
    dz = cell_size[2] / nz if nz > 0 else 1.0

    nx1, ny1, nz1 = nx + 1, ny + 1, nz + 1
    total_grid = nx1 * ny1 * nz1

    id_to_idx = mesh_data['vertices']
    positions = {}
    for nid in id_to_idx:
        idx = nid - 1
        i = idx % nx1
        j = (idx // nx1) % ny1
        k = idx // (nx1 * ny1)
        x = -cell_size[0] / 2 + i * dx
        y = -cell_size[1] / 2 + j * dy
        z = -cell_size[2] / 2 + k * dz
        positions[nid] = np.array([x, y, z]) * length_scale

    return positions


def write_obj(mesh_data, output_path, vertex_positions=None):
    id_to_idx = mesh_data['vertices']
    inv_map = {v: k for k, v in id_to_idx.items()}

    with open(output_path, 'w') as f:
        f.write(f"# Generated from RecurDyn RMD GGEOM\n")
        f.write(f"# patches: {mesh_data['num_patches']}, nodes: {mesh_data['num_nodes']}\n")
        f.write(f"# extracted triangles: {mesh_data['num_triangles']}\n")

        for local_idx in range(1, len(id_to_idx) + 1):
            nid = inv_map[local_idx]
            if vertex_positions is not None and nid in vertex_positions:
                p = vertex_positions[nid]
                f.write(f"v {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
            else:
                f.write(f"v 0 0 0\n")

        for t in mesh_data['triangles']:
            f.write(f"f {t[0]} {t[1]} {t[2]}\n")

    return output_path


def extract_to_obj(entities, ggeom_id, output_dir, label='surface', length_scale=1.0):
    mesh = extract_mesh_from_rmd(entities, ggeom_id)
    if mesh is None:
        print(f"  No mesh data found for GGEOM {ggeom_id}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    obj_path = os.path.join(output_dir, f"{label}_ggeom{ggeom_id}.obj")

    positions = assign_grid_positions(mesh, length_scale)
    write_obj(mesh, obj_path, positions)

    print(f"  GGEOM {ggeom_id}: {mesh['num_triangles']} triangles, "
          f"{mesh['num_unique_vertices']} vertices -> {obj_path}")
    return obj_path


def extract_all_surfaces(entities, output_dir, length_scale=1.0):
    paths = {}
    for ent in entities:
        if ent.kind == 'GGEOM':
            npatch = ent.get_int('NO_PATCH')
            label = ent.get_str('NAME', f'surface_{ent.id}').replace('##GSURFACE##_', '')
            path = extract_to_obj(entities, ent.id, output_dir, label, length_scale)
            if path:
                paths[ent.id] = path
    return paths
