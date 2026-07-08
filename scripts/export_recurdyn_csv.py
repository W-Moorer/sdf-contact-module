"""
Export RecurDyn .rplt results to CSV via COM API.

Usage: python scripts/export_recurdyn_csv.py

Requirements:
  - RecurDyn 2023 installed
  - pywin32 (pip install pywin32)
"""

import os, sys, csv
import win32com.client, pythoncom, pywintypes

RD_CLSID = '{FEAC72DD-1D1B-47D0-8958-5BF3E69F0A93}'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'dir': os.path.join('..', 'models', '圆环-立方体-对心碰撞'),
        'rdyn': 'Ring - Cube - Centric Collision.rdyn',
        'result_subdir': 'Ring - Cube - Centric Collision_01',
        'signals': ['Bodies/Body1/Pos_TZ',
                    'Bodies/Body1/Pos_TX', 'Bodies/Body1/Pos_TY',
                    'Bodies/Body1/Vel_TZ',
                    'Bodies/Body2/Pos_TZ',
                    'Contacts/GeoSurContact1/Fn_Total'],
    },
    {
        'label': 'gear_transmission_01',
        'dir': os.path.join('..', 'models', '齿轮-齿轮-一级转动传动'),
        'rdyn': 'First-stage gear transmission.rdyn',
        'result_subdir': 'First-stage gear transmission_01',
        'signals': ['Bodies/Verbindungsrad_Rechts6/W_RZ',
                    'Bodies/Verbindungsrad_Links6/W_RZ',
                    'Joints/RevJoint1/Displacement_Z',
                    'Joints/RevJoint2/Displacement_Z',
                    'Contacts/GeoSurContact1/Fn_Total',
                    'Contacts/GeoSurContact1/Ff_Total'],
    },
    {
        'label': 'hollow_cyl_friction',
        'dir': os.path.join('..', 'models', '空心圆柱-圆柱-摩擦接触'),
        'rdyn': 'Hollow cylinder - Cylinder - Frictional contact.rdyn',
        'result_subdir': 'Hollow cylinder - Cylinder - Frictional contact_01',
        'signals': ['Bodies/Body1/W_RZ',
                    'Bodies/Body2/W_RZ',
                    'Joints/Cylindrical1/Torque_Z',
                    'Contacts/GeoSurContact1/Fn_Total',
                    'Contacts/GeoSurContact1/Ff_Total'],
    },
]


def find_rplt(scenario):
    base = os.path.abspath(os.path.join(SCRIPT_DIR, scenario['dir']))
    rdyn_path = os.path.join(base, scenario['rdyn'])
    result_dir = os.path.join(base, scenario['result_subdir'])
    rplt_files = [f for f in os.listdir(result_dir) if f.endswith('.rplt')] if os.path.isdir(result_dir) else []
    rplt_path = os.path.join(result_dir, rplt_files[0]) if rplt_files else None
    return rdyn_path, rplt_path


def init_com():
    try:
        return win32com.client.Dispatch("RecurDyn.Application")
    except:
        clsid = pywintypes.IID(RD_CLSID)
        obj = pythoncom.CoCreateInstance(clsid, None,
                                         pythoncom.CLSCTX_LOCAL_SERVER,
                                         pythoncom.IID_IDispatch)
        return win32com.client.Dispatch(obj).RecurDynApplication


def export_scenario(rdApp, scenario):
    label = scenario['label']
    rdyn_path, rplt_path = find_rplt(scenario)

    if not os.path.exists(rdyn_path):
        print(f"  [{label}] SKIP: no .rdyn")
        return False
    if not rplt_path or not os.path.exists(rplt_path):
        print(f"  [{label}] SKIP: no .rplt")
        return False

    print(f"\n  [{label}]")
    print(f"    RPLT: {rplt_path}")

    out_dir = os.path.join(OUTPUT_DIR, label)
    os.makedirs(out_dir, exist_ok=True)

    # Copy RPLT to temp to avoid path issues
    temp_dir = os.path.join(os.environ['TEMP'], 'rd_csv_export', label)
    os.makedirs(temp_dir, exist_ok=True)
    import shutil
    temp_rplt = os.path.join(temp_dir, 'results.rplt')
    shutil.copy2(rplt_path, temp_rplt)

    # If model not open, open it to get a document context
    doc = rdApp.OpenModelDocument(rdyn_path)

    try:
        # Create empty plot document
        plot_doc = doc.CreatePlotDocument(1)  # PlotDocType.Empty = 1
        print(f"    PlotDocument created")

        # Import RPLT
        idx = plot_doc.ImportPlotDataFile(temp_rplt)
        print(f"    RPLT imported, index={idx}")

        # Export each signal to CSV
        for signal_path in scenario['signals']:
            try:
                data = plot_doc.GetPlotDataFromIndex(idx, signal_path)
                if data:
                    time_data = plot_doc.GetPlotDataFromIndex(idx, 'Time')
                    csv_name = signal_path.replace('/', '_') + '.csv'
                    csv_path = os.path.join(out_dir, csv_name)
                    with open(csv_path, 'w', newline='') as f:
                        w = csv.writer(f)
                        w.writerow(['time', 'value'])
                        for t, v in zip(time_data, data):
                            w.writerow([t, v])
                    print(f"    {csv_name}: {len(data)} pts")
                else:
                    print(f"    {signal_path}: empty")
            except Exception as ex:
                print(f"    {signal_path}: {str(ex)[:60]}")

        print(f"    Done -> {out_dir}")
        return True

    except Exception as ex:
        import traceback
        traceback.print_exc()
        return False
    finally:
        try:
            rdApp.CloseAllDocument()
        except:
            pass


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    rdApp = init_com()
    print(f"RecurDyn Application ready")

    for scenario in SCENARIOS:
        export_scenario(rdApp, scenario)

    print(f"\nAll exports: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
