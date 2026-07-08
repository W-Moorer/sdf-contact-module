"""
Export RecurDyn results to CSV by:
  A) Copying existing .out/.req files (already have time-stepping data)
  B) Writing a RecurDyn macro script for in-application export
  C) Trying COM API for direct export

Usage:
    python scripts/export_all_csv.py
"""

import os, sys, shutil, glob

SCRIPT_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'result_dir': os.path.join('..', 'models', '圆环-立方体-对心碰撞',
                                   'Ring - Cube - Centric Collision_01'),
        'files': ['out', 'req', 'rad', 'rdp'],
    },
    {
        'label': 'gear_transmission_01',
        'result_dir': os.path.join('..', 'models', '齿轮-齿轮-一级转动传动',
                                   'First-stage gear transmission_01'),
        'files': ['out', 'req', 'rad', 'rdp'],
    },
    {
        'label': 'gear_transmission_02',
        'result_dir': os.path.join('..', 'models', '齿轮-齿轮-一级转动传动',
                                   'First-stage gear transmission_02'),
        'files': ['out', 'req', 'rad', 'rdp'],
    },
    {
        'label': 'gear_transmission_03',
        'result_dir': os.path.join('..', 'models', '齿轮-齿轮-一级转动传动',
                                   'First-stage gear transmission_03'),
        'files': ['out', 'req', 'rad', 'rdp'],
    },
    {
        'label': 'hollow_cyl_friction',
        'result_dir': os.path.join('..', 'models', '空心圆柱-圆柱-摩擦接触',
                                   'Hollow cylinder - Cylinder - Frictional contact_01'),
        'files': ['out', 'req', 'rad', 'rdp'],
    },
]


def copy_existing_results():
    """Copy all existing ASCII result files to reference_csv/"""
    for cfg in SCENARIOS:
        src = os.path.abspath(os.path.join(SCRIPT_DIR, cfg['result_dir']))
        dst = os.path.join(OUTPUT_DIR, cfg['label'])
        os.makedirs(dst, exist_ok=True)
        
        if not os.path.exists(src):
            print(f"  SKIP {cfg['label']}: {src} not found")
            continue
        
        count = 0
        for fname in os.listdir(src):
            fpath = os.path.join(src, fname)
            if os.path.isfile(fpath):
                ext = fname.split('.')[-1].lower() if '.' in fname else ''
                if ext in cfg['files']:
                    shutil.copy2(fpath, os.path.join(dst, fname))
                    count += 1
        print(f"  {cfg['label']}: copied {count} files")

    # Also check .rss files
    rss_gear = os.path.join('..', 'models', '齿轮-齿轮-一级转动传动',
                            'run_dynamic_1s_1000step.rss')
    rss_path = os.path.abspath(os.path.join(SCRIPT_DIR, rss_gear))
    if os.path.exists(rss_path):
        shutil.copy2(rss_path, os.path.join(OUTPUT_DIR, 'gear_transmission_01', 'run.rss'))
        print(f"  Copied .rss file")


def write_recurdyn_macro():
    """Write a Python macro to run inside RecurDyn for CSV export."""
    macro_path = os.path.join(SCRIPT_DIR, 'export_from_recurdyn.py')
    with open(macro_path, 'w') as f:
        f.write(r'''"""
RecurDyn Python Macro: Export simulation results to CSV.

HOW TO USE:
1. Open RecurDyn 2023
2. Tools -> Python -> Run Script
3. Select this file
4. CSV files will be saved to reference_csv/<label>/

This macro opens each .rdyn project, runs the solver,
and exports key signals to CSV.
"""

import os, csv, sys

# ---- Configuration ----
BASE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(BASE)
OUTPUT = os.path.join(PARENT, 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'rdyn': os.path.join(PARENT, 'models', '圆环-立方体-对心碰撞',
                             'Ring - Cube - Centric Collision.rdyn'),
        'export': [
            ('Body1_Position', 'body_ring_pos'),
            ('Body1_Velocity', 'body_ring_vel'),
            ('Body2_Position', 'body_cube_vel'),
            ('GeoSurContact1_Force', 'contact_force'),
        ],
    },
    {
        'label': 'gear_transmission_01',
        'rdyn': os.path.join(PARENT, 'models', '齿轮-齿轮-一级转动传动',
                             'First-stage gear transmission.rdyn'),
        'rss': os.path.join(PARENT, 'models', '齿轮-齿轮-一级转动传动',
                            'run_dynamic_1s_1000step.rss'),
        'export': [
            ('RevJoint1_Displacement', 'joint_driving_theta'),
            ('RevJoint2_Displacement', 'joint_driven_theta'),
            ('Verbindungsrad_Rechts6_AngularVelocity', 'body_driving_omega'),
            ('Verbindungsrad_Links6_AngularVelocity', 'body_driven_omega'),
        ],
    },
    {
        'label': 'hollow_cyl_friction',
        'rdyn': os.path.join(PARENT, 'models', '空心圆柱-圆柱-摩擦接触',
                             'Hollow cylinder - Cylinder - Frictional contact.rdyn'),
        'export': [
            ('Body1_Position', 'body_hollow_pos'),
            ('Body1_AngularVelocity', 'body_hollow_omega'),
            ('Cylindrical1_Torque', 'joint_friction_torque'),
            ('GeoSurContact1_Force', 'contact_force'),
        ],
    },
]


def export_signal(model, signal_name, output_path, nsteps):
    """Export a signal from the model results to CSV."""
    try:
        result = model.Results(signal_name)
    except:
        print(f"    Signal not found: {signal_name}")
        return False
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['time', 'value'])
        for i in range(nsteps + 1):
            t = result.GetTime(i)
            v = result.GetValue(i)
            writer.writerow([t, v])
    print(f"    Exported: {signal_name} -> {os.path.basename(output_path)}")
    return True


def main():
    print("=" * 60)
    print("  RecurDyn CSV Export Macro")
    print("=" * 60)
    
    rdApp = RecurDyn.Application
    rdApp.Visible = False
    
    for cfg in SCENARIOS:
        label = cfg['label']
        rdyn = cfg['rdyn']
        
        print(f"\n  [{label}]")
        print(f"  Opening: {rdyn}")
        
        if not os.path.exists(rdyn):
            print(f"  SKIP: file not found")
            continue
        
        doc = rdApp.OpenModelDocument(rdyn)
        model = doc.Model
        
        rss = cfg.get('rss')
        if rss and os.path.exists(rss):
            try:
                doc.LoadSimulationSetting(rss)
                print(f"  Loaded .rss settings")
            except:
                pass
        
        print(f"  Solving...")
        doc.Solve()
        print(f"  Solved")
        
        out_dir = os.path.join(OUTPUT, label)
        for sig_name, file_label in cfg['export']:
            csv_path = os.path.join(out_dir, f"{file_label}.csv")
            export_signal(model, sig_name, csv_path, 1000)
        
        rdApp.CloseAllDocument()
    
    print(f"\n  Done. CSV files in: {OUTPUT}")


if __name__ == '__main__':
    main()
''')
    print(f"  Macro written: {macro_path}")
    return macro_path


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("Step A: Copy existing result files...")
    copy_existing_results()
    
    print("\nStep B: Write RecurDyn macro...")
    macro = write_recurdyn_macro()
    
    print(f"\nTo export CSV data from RecurDyn GUI:")
    print(f"  1. Open RecurDyn 2023")
    print(f"  2. Tools -> Python -> Run Script")
    print(f"  3. Select: {macro}")
    print(f"  4. CSV files saved to: {os.path.abspath(OUTPUT_DIR)}")


if __name__ == '__main__':
    main()
