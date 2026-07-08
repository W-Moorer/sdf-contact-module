# scripts — Utility Scripts

| Script | Description |
|--------|-------------|
| `gen_meshes.py` | Generate OBJ meshes (cube, hollow cylinder variants) for SDF contact validation |
| `export_all_signals.py` | Batch-export RecurDyn `.rplt` → CSV via COM API + PlotDocument |
| `export_all_data.py` | Export all 30 signals per scenario from RPLT |
| `compare_and_plot.py` | Three-way overlay plots: RecurDyn vs Framework vs Analytical |
| `export_from_recurdyn.py` | RecurDyn internal macro for CSV export (run from RecurDyn GUI) |

## COM API Export

Requires RecurDyn 2023 installed and COM DLL registered:
```
regsvr32 "C:\Program Files\FunctionBay, Inc\RecurDyn 2023\Bin\RecurDynCOM.dll"
```
