import os


def extract_mesh_from_rmd(entities, ggeom_id):
    ggeoms = [e for e in entities if e.kind == 'GGEOM' and e.id == ggeom_id]
    if not ggeoms:
        return None
    ent = ggeoms[0]

    npatch = ent.get_int('NO_PATCH')
    nnodes = ent.get_int('NO_NODE')
    ncell = ent.get_floats('NCELL')
    cell_size = ent.get_floats('CELL_SIZE')

    node_positions = {}
    triangles = []
    mode = None

    for line in ent.raw_lines:
        s = line.strip()
        if not s or s.startswith('!'):
            continue

        if ', NODES' in s.upper():
            mode = 'nodes'
            continue
        if ', PATCHES' in s.upper():
            mode = 'patches'
            continue

        if mode == 'nodes':
            parts = [p.strip() for p in s.lstrip(',').split(',') if p.strip()]
            if len(parts) == 3:
                try:
                    x, y, z = map(float, parts)
                    node_positions[len(node_positions) + 1] = (x, y, z)
                except ValueError:
                    pass

        elif mode == 'patches' and s.startswith('3,'):
            parts = s.split(',')
            if len(parts) >= 4:
                try:
                    n1, n2, n3 = int(parts[1]), int(parts[2]), int(parts[3])
                    triangles.append((n1, n2, n3))
                except ValueError:
                    pass

    if len(node_positions) == 0 or len(triangles) == 0:
        return None

    return {
        'num_patches': npatch,
        'num_nodes': nnodes,
        'num_triangles': len(triangles),
        'num_unique_vertices': len(node_positions),
        'node_positions': node_positions,
        'triangles': triangles,
        'ncell': ncell,
        'cell_size': cell_size,
    }


def write_obj(mesh_data, output_path, length_scale=1.0):
    with open(output_path, 'w') as f:
        f.write(f"# Generated from RecurDyn RMD GGEOM\n")
        f.write(f"# patches: {mesh_data['num_patches']}, nodes: {mesh_data['num_nodes']}\n")
        f.write(f"# extracted triangles: {len(mesh_data['triangles'])}\n\n")

        for nid in sorted(mesh_data['node_positions']):
            x, y, z = mesh_data['node_positions'][nid]
            f.write(f"v {x*length_scale:.10f} {y*length_scale:.10f} {z*length_scale:.10f}\n")

        f.write('\n')
        for n1, n2, n3 in mesh_data['triangles']:
            f.write(f"f {n1} {n2} {n3}\n")

    return output_path


def extract_to_obj(entities, ggeom_id, output_dir, label='surface', length_scale=1.0):
    mesh = extract_mesh_from_rmd(entities, ggeom_id)
    if mesh is None:
        print(f"  No mesh data found for GGEOM {ggeom_id}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    obj_path = os.path.join(output_dir, f"{label}_ggeom{ggeom_id}.obj")
    write_obj(mesh, obj_path, length_scale)

    print(f"  GGEOM {ggeom_id}: {mesh['num_triangles']} tris, "
          f"{mesh['num_unique_vertices']} verts -> {os.path.basename(obj_path)}")
    return obj_path


def extract_all_surfaces(entities, output_dir, length_scale=1.0):
    paths = {}
    for ent in entities:
        if ent.kind == 'GGEOM':
            label = ent.get_str('NAME', f'surface_{ent.id}').replace('##GSURFACE##_', '')
            path = extract_to_obj(entities, ent.id, output_dir, label, length_scale)
            if path:
                paths[ent.id] = path
    return paths
