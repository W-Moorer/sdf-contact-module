"""Export ALL signals from RPLT using GetPlottableNameList."""
import os, csv, shutil, win32com.client, pythoncom, pywintypes

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPTS_DIR, '..', 'reference_csv')
RD_CLSID = '{FEAC72DD-1D1B-47D0-8958-5BF3E69F0A93}'

SCENARIOS = [
    ('hollow_cyl_friction', '..\\models\\空心圆柱-圆柱-摩擦接触',
     'Hollow cylinder - Cylinder - Frictional contact.rdyn',
     'Hollow cylinder - Cylinder - Frictional contact_01'),
    ('ring_cube_collision', '..\\models\\圆环-立方体-对心碰撞',
     'Ring - Cube - Centric Collision.rdyn',
     'Ring - Cube - Centric Collision_01'),
    ('gear_transmission_01', '..\\models\\齿轮-齿轮-一级转动传动',
     'First-stage gear transmission.rdyn',
     'First-stage gear transmission_01'),
    ('gear_transmission_02', '..\\models\\齿轮-齿轮-一级转动传动',
     'First-stage gear transmission.rdyn',
     'First-stage gear transmission_02'),
    ('gear_transmission_03', '..\\models\\齿轮-齿轮-一级转动传动',
     'First-stage gear transmission.rdyn',
     'First-stage gear transmission_03'),
]

clsid = pywintypes.IID(RD_CLSID)
obj = pythoncom.CoCreateInstance(clsid, None, pythoncom.CLSCTX_LOCAL_SERVER, pythoncom.IID_IDispatch)
rdApp = win32com.client.Dispatch(obj).RecurDynApplication

for label, rdyn_dir, rdyn_name, result_subdir in SCENARIOS:
    base = os.path.abspath(os.path.join(SCRIPTS_DIR, rdyn_dir))
    rdyn = os.path.join(base, rdyn_name)
    result_dir = os.path.join(base, result_subdir)
    rplt_list = [f for f in os.listdir(result_dir) if f.endswith('.rplt')]
    if not rplt_list:
        print(f'{label}: no RPLT')
        continue

    rplt = os.path.join(result_dir, rplt_list[0])
    temp = os.path.join(os.environ['TEMP'], 'rd_export_full', label)
    os.makedirs(temp, exist_ok=True)
    shutil.copy2(rplt, os.path.join(temp, 'r.rplt'))

    doc = rdApp.OpenModelDocument(rdyn)
    pd = doc.CreatePlotDocument(1)
    idx = pd.ImportPlotDataFile(os.path.join(temp, 'r.rplt'))

    names = pd.GetPlottableNameList(idx)
    time_data = pd.GetPlotDataFromIndex(idx, 'TIME')

    out_dir = os.path.join(OUTPUT_DIR, label)
    os.makedirs(out_dir, exist_ok=True)

    count = 0
    for sig in names:
        if sig == 'TIME':
            continue
        try:
            data = pd.GetPlotDataFromIndex(idx, sig)
            if data and len(data) > 0:
                safe = sig.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
                csv_path = os.path.join(out_dir, f'{safe}.csv')
                with open(csv_path, 'w', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(['time', 'value'])
                    for t, v in zip(time_data, data):
                        w.writerow([t, v])
                count += 1
        except:
            pass

    print(f'[{label}] {count}/{len(names)} signals exported')
    rdApp.CloseAllDocument()
