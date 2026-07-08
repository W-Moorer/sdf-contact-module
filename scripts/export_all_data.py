"""
Final comprehensive export: ALL signals from ALL 5 RecurDyn scenarios.
Then run comparison against our framework + analytical solutions.
"""
import os, csv, shutil, sys
import win32com.client, pythoncom, pywintypes
import numpy as np

RD_CLSID = '{FEAC72DD-1D1B-47D0-8958-5BF3E69F0A93}'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'path': '..\\models\\圆环-立方体-对心碰撞',
        'rdyn': 'Ring - Cube - Centric Collision.rdyn',
        'sub': 'Ring - Cube - Centric Collision_01',
        'entities': {'Bodies': ['Body1', 'Body2']},
    },
    {
        'label': 'gear_transmission_01',
        'path': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'sub': 'First-stage gear transmission_01',
        'entities': {'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6']},
    },
    {
        'label': 'gear_transmission_02',
        'path': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'sub': 'First-stage gear transmission_02',
        'entities': {'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6']},
    },
    {
        'label': 'gear_transmission_03',
        'path': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'sub': 'First-stage gear transmission_03',
        'entities': {'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6']},
    },
    {
        'label': 'hollow_cyl_friction',
        'path': '..\\models\\空心圆柱-圆柱-摩擦接触',
        'rdyn': 'Hollow cylinder - Cylinder - Frictional contact.rdyn',
        'sub': 'Hollow cylinder - Cylinder - Frictional contact_01',
        'entities': {'Bodies': ['Body1', 'Body2']},
    },
]

COMPONENTS = ['Pos_TX', 'Pos_TY', 'Pos_TZ',
              'Vel_TX', 'Vel_TY', 'Vel_TZ',
              'Vel_RX', 'Vel_RY', 'Vel_RZ',
              'Acc_TX', 'Acc_TY', 'Acc_TZ',
              'Acc_RX', 'Acc_RY', 'Acc_RZ']


def export_all():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    clsid = pywintypes.IID(RD_CLSID)
    obj = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER, pythoncom.IID_IDispatch)
    rdApp = win32com.client.Dispatch(obj).RecurDynApplication

    for sc in SCENARIOS:
        label = sc['label']
        base = os.path.abspath(os.path.join(SCRIPT_DIR, sc['path']))
        rdyn_path = os.path.join(base, sc['rdyn'])
        result_dir = os.path.join(base, sc['sub'])
        rplt_list = [f for f in os.listdir(result_dir) if f.endswith('.rplt')]
        if not rplt_list:
            print(f'[{label}] no RPLT')
            continue

        temp = os.path.join(os.environ['TEMP'], 'rd_final', label)
        os.makedirs(temp, exist_ok=True)
        shutil.copy2(os.path.join(result_dir, rplt_list[0]), os.path.join(temp, 'r.rplt'))

        try:
            doc = rdApp.OpenModelDocument(rdyn_path)
            pd = doc.CreatePlotDocument(1)
            idx = pd.ImportPlotDataFile(os.path.join(temp, 'r.rplt'))
            time_data = pd.GetPlotDataFromIndex(idx, 'Time')
            dt = time_data[1] - time_data[0] if len(time_data) > 1 else 0.001
            print(f'[{label}] dt={dt:.6f}s, steps={len(time_data)}')

            out_dir = os.path.join(OUTPUT_DIR, label)
            os.makedirs(out_dir, exist_ok=True)

            count = 0
            for cat, entities in sc['entities'].items():
                for ent in entities:
                    for comp in COMPONENTS:
                        sig = f'{cat}/{ent}/{comp}'
                        try:
                            data = pd.GetPlotDataFromIndex(idx, sig)
                            if data and len(data) > 0:
                                name = sig.replace('/', '_')
                                with open(os.path.join(out_dir, f'{name}.csv'), 'w', newline='') as f:
                                    w = csv.writer(f)
                                    w.writerow(['time', 'value'])
                                    for t, v in zip(time_data, data):
                                        w.writerow([t, v])
                                count += 1
                        except:
                            pass

            # Also copy original files
            for ext in ['.req', '.out', '.rss']:
                for f in os.listdir(result_dir):
                    if f.endswith(ext):
                        shutil.copy2(os.path.join(result_dir, f), os.path.join(out_dir, f))

            print(f'  -> {count} signals')
            rdApp.CloseAllDocument()
        except Exception as ex:
            print(f'[{label}] ERROR: {ex}')
            try: rdApp.CloseAllDocument()
            except: pass


if __name__ == '__main__':
    export_all()
    print(f'\nAll data exported to {OUTPUT_DIR}')
