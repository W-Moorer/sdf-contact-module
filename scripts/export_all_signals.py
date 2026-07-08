"""
Export all available signals from RecurDyn .rplt files to CSV.
"""
import os, csv, shutil
import win32com.client, pythoncom, pywintypes

RD_CLSID = '{FEAC72DD-1D1B-47D0-8958-5BF3E69F0A93}'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'reference_csv')

SCENARIOS = [
    {
        'label': 'ring_cube_collision',
        'rdyn_dir': '..\\models\\圆环-立方体-对心碰撞',
        'rdyn': 'Ring - Cube - Centric Collision.rdyn',
        'result_subdir': 'Ring - Cube - Centric Collision_01',
        'entities': {
            'Bodies': ['Body1', 'Body2'],
        },
    },
    {
        'label': 'gear_transmission_01',
        'rdyn_dir': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'result_subdir': 'First-stage gear transmission_01',
        'entities': {
            'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6'],
        },
    },
    {
        'label': 'gear_transmission_02',
        'rdyn_dir': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'result_subdir': 'First-stage gear transmission_02',
        'entities': {
            'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6'],
        },
    },
    {
        'label': 'gear_transmission_03',
        'rdyn_dir': '..\\models\\齿轮-齿轮-一级转动传动',
        'rdyn': 'First-stage gear transmission.rdyn',
        'result_subdir': 'First-stage gear transmission_03',
        'entities': {
            'Bodies': ['Verbindungsrad_Rechts6', 'Verbindungsrad_Links6'],
        },
    },
    {
        'label': 'hollow_cyl_friction',
        'rdyn_dir': '..\\models\\空心圆柱-圆柱-摩擦接触',
        'rdyn': 'Hollow cylinder - Cylinder - Frictional contact.rdyn',
        'result_subdir': 'Hollow cylinder - Cylinder - Frictional contact_01',
        'entities': {
            'Bodies': ['Body1', 'Body2'],
        },
    },
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    clsid = pywintypes.IID(RD_CLSID)
    obj = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER, pythoncom.IID_IDispatch)
    rdApp = win32com.client.Dispatch(obj).RecurDynApplication
    
    components = {
        'Pos_TX': 'Position_X', 'Pos_TY': 'Position_Y', 'Pos_TZ': 'Position_Z',
        'Vel_TX': 'Velocity_X', 'Vel_TY': 'Velocity_Y', 'Vel_TZ': 'Velocity_Z',
        'Acc_TX': 'Accel_X', 'Acc_TY': 'Accel_Y', 'Acc_TZ': 'Accel_Z',
    }
    
    for sc in SCENARIOS:
        label = sc['label']
        rdyn_base = os.path.abspath(os.path.join(SCRIPT_DIR, sc['rdyn_dir']))
        rdyn_path = os.path.join(rdyn_base, sc['rdyn'])
        result_dir = os.path.join(rdyn_base, sc['result_subdir'])
        rplt_files = [f for f in os.listdir(result_dir) if f.endswith('.rplt')]
        
        if not rplt_files:
            print(f'[{label}] no .rplt, SKIP')
            continue
        
        rplt_path = os.path.join(result_dir, rplt_files[0])
        
        # Copy to temp
        temp = os.path.join(os.environ['TEMP'], 'rd_export_all', label)
        os.makedirs(temp, exist_ok=True)
        shutil.copy2(rplt_path, os.path.join(temp, 'results.rplt'))
        
        try:
            doc = rdApp.OpenModelDocument(rdyn_path)
            plot_doc = doc.CreatePlotDocument(1)
            idx = plot_doc.ImportPlotDataFile(os.path.join(temp, 'results.rplt'))
            
            out_dir = os.path.join(OUTPUT_DIR, label)
            os.makedirs(out_dir, exist_ok=True)
            
            # Get time vector once
            time_data = plot_doc.GetPlotDataFromIndex(idx, 'Time')
            
            count = 0
            for cat, entities in sc['entities'].items():
                for ent in entities:
                    for comp_rplt, comp_name in components.items():
                        sig = f'{cat}/{ent}/{comp_rplt}'
                        try:
                            data = plot_doc.GetPlotDataFromIndex(idx, sig)
                            if data and len(data) > 0:
                                csv_name = f'{cat}_{ent}_{comp_name}.csv'
                                csv_path = os.path.join(out_dir, csv_name)
                                with open(csv_path, 'w', newline='') as f:
                                    w = csv.writer(f)
                                    w.writerow(['time', 'value'])
                                    for t, v in zip(time_data, data):
                                        w.writerow([t, v])
                                count += 1
                        except:
                            pass
            
            # Also copy existing .req, .out files
            for ext in ['.req', '.out', '.rss', '.rdp']:
                for f in os.listdir(result_dir):
                    if f.endswith(ext):
                        shutil.copy2(os.path.join(result_dir, f), os.path.join(out_dir, f))
            
            print(f'[{label}] {count} signals -> {out_dir}')
            rdApp.CloseAllDocument()
        except Exception as ex:
            print(f'[{label}] ERROR: {ex}')
            try: rdApp.CloseAllDocument()
            except: pass
    
    print(f'\nAll done. Files in: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
