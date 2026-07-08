"""
RecurDyn Python Macro: Export all 5 scenarios to CSV.

HOW TO USE:
  1. Open RecurDyn 2023
  2. Tools -> Python -> Run Script
  3. Select this file
  4. CSV files saved to ../reference_csv/<label>/
"""

import os, sys, csv, win32com.client

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT = os.path.join(ROOT, 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'rdyn': os.path.join(ROOT, 'models', '圆环-立方体-对心碰撞',
                             'Ring - Cube - Centric Collision.rdyn'),
    },
    {
        'label': 'gear_transmission_01',
        'rdyn': os.path.join(ROOT, 'models', '齿轮-齿轮-一级转动传动',
                             'First-stage gear transmission.rdyn'),
    },
    {
        'label': 'gear_transmission_02',
        'rdyn': os.path.join(ROOT, 'models', '齿轮-齿轮-一级转动传动',
                             'First-stage gear transmission.rdyn'),
    },
    {
        'label': 'gear_transmission_03',
        'rdyn': os.path.join(ROOT, 'models', '齿轮-齿轮-一级转动传动',
                             'First-stage gear transmission.rdyn'),
    },
    {
        'label': 'hollow_cyl_friction',
        'rdyn': os.path.join(ROOT, 'models', '空心圆柱-圆柱-摩擦接触',
                             'Hollow cylinder - Cylinder - Frictional contact.rdyn'),
    },
]


def export_model(rdApp, rdyn_path, label):
    print(f"\n--- {label} ---")
    print(f"Opening: {rdyn_path}")
    doc = rdApp.OpenModelDocument(rdyn_path)
    model = doc.Model
    print(f"Model: {model.Name}")
    
    # Run solver
    print("Solving...")
    doc.AnalysisStart()
    print("Solved")
    
    out_dir = os.path.join(OUTPUT, label)
    os.makedirs(out_dir, exist_ok=True)
    
    # Export all body results
    for body_idx in range(1, 10):
        try:
            body = model.Bodies.Item(body_idx)
            bname = body.Name.replace(' ', '_')
            
            for comp, axis in [('Position', 'X'), ('Position', 'Y'), ('Position', 'Z'),
                               ('Velocity', 'X'), ('Velocity', 'Y'), ('Velocity', 'Z'),
                               ('AngularVelocity', 'X'), ('AngularVelocity', 'Y'), ('AngularVelocity', 'Z')]:
                try:
                    data = doc.GetPlotData(body.Name, comp, axis)
                    if data:
                        times = doc.GetPlotDataTime(body.Name, comp, axis)
                        csv_path = os.path.join(out_dir, f"body_{bname}_{comp}_{axis}.csv")
                        with open(csv_path, 'w', newline='') as f:
                            w = csv.writer(f)
                            w.writerow(['time', 'value'])
                            for t, v in zip(times, data):
                                w.writerow([t, v])
                        print(f"  body_{bname}_{comp}_{axis}: {len(data)} pts")
                except:
                    pass
        except:
            break
    
    # Export joint data
    for joint_idx in range(1, 10):
        try:
            joint = model.Joints.Item(joint_idx)
            jname = joint.Name.replace(' ', '_')
            
            for comp in ['Displacement', 'Velocity', 'Force', 'Torque']:
                for axis in ['X', 'Y', 'Z', 'M']:
                    try:
                        data = doc.GetPlotData(joint.Name, comp, axis)
                        if data:
                            times = doc.GetPlotDataTime(joint.Name, comp, axis)
                            csv_path = os.path.join(out_dir, f"joint_{jname}_{comp}_{axis}.csv")
                            with open(csv_path, 'w', newline='') as f:
                                w = csv.writer(f)
                                w.writerow(['time', 'value'])
                                for t, v in zip(times, data):
                                    w.writerow([t, v])
                            print(f"  joint_{jname}_{comp}_{axis}: {len(data)} pts")
                    except:
                        pass
        except:
            break
    
    # Export contact data
    for ct_idx in range(1, 10):
        try:
            ct = model.GeoSurfaceContacts.Item(ct_idx)
            ctname = ct.Name.replace(' ', '_')
            
            for comp in ['Force', 'Penetration', 'FrictionForce']:
                for axis in ['X', 'Y', 'Z', 'M']:
                    try:
                        data = doc.GetPlotData(ct.Name, comp, axis)
                        if data:
                            times = doc.GetPlotDataTime(ct.Name, comp, axis)
                            csv_path = os.path.join(out_dir, f"contact_{ctname}_{comp}_{axis}.csv")
                            with open(csv_path, 'w', newline='') as f:
                                w = csv.writer(f)
                                w.writerow(['time', 'value'])
                                for t, v in zip(times, data):
                                    w.writerow([t, v])
                            print(f"  contact_{ctname}_{comp}_{axis}: {len(data)} pts")
                    except:
                        pass
        except:
            break
    
    rdApp.CloseAllDocument()
    print(f"Done: {label}")


def main():
    rdApp = win32com.client.Dispatch("RecurDyn.Application")
    
    for cfg in SCENARIOS:
        rdyn = cfg['rdyn']
        if os.path.exists(rdyn):
            export_model(rdApp, rdyn, cfg['label'])
        else:
            print(f"\nSKIP {cfg['label']}: {rdyn} not found")
    
    print(f"\nAll exports done. CSV files in: {OUTPUT}")

if __name__ == '__main__':
    main()
