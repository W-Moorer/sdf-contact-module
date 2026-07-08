import os
import numpy as np


class ValidationCase:
    def __init__(self, name, rmd_path, T, dt,
                 reference_signals=None, tolerances=None):
        self.name = name
        self.rmd_path = rmd_path
        self.T = T
        self.dt = dt
        self.reference_signals = reference_signals or {}
        self.tolerances = tolerances or {'smoke': 0.25, 'strict': 0.1}


class RecurDynValidator:
    def __init__(self):
        self.results = []

    def run_case(self, case, model, history):
        from .signal_metrics import compute_error_metrics, classify_error

        result = {
            'case_name': case.name,
            'T': case.T,
            'dt': case.dt,
            'steps': len(history),
            'signals': {},
        }

        for signal_name, ref in case.reference_signals.items():
            our_signal = self._extract_signal(history, signal_name)
            if our_signal is None:
                continue

            ref_values = np.array(ref['values'])
            ref_time = np.array(ref['time'])

            our_time = np.array([r.get('t', 0) for r in history])

            metrics = compute_error_metrics(ref_values, our_signal, ref_time, our_time)
            status, level = classify_error(metrics,
                                           case.tolerances.get('smoke', 0.25),
                                           case.tolerances.get('strict', 0.1))

            result['signals'][signal_name] = {
                'metrics': metrics,
                'status': status,
                'level': level,
                'ref_last': float(ref_values[-1]),
                'our_last': float(our_signal[-1]),
            }

        self.results.append(result)
        return result

    def _extract_signal(self, history, signal_name):
        parts = signal_name.split('.')
        if len(parts) == 2 and parts[0] == 'joint':
            key = f'joint_{parts[1]}_theta'
            return np.array([r.get(key, np.nan) for r in history])
        if len(parts) == 2 and parts[0] == 'body':
            if parts[1].endswith('_omega_z'):
                body_name = parts[1].replace('_omega_z', '')
                key_r = f'body_{body_name}_r'
                omegas = []
                for r in history:
                    omega = r.get(key_r.replace('_r', '_v'), None)
                    omegas.append(np.nan)
                return np.array(omegas)
            key_r = f'body_{parts[1]}_r'
            return np.array([r.get(key_r, np.zeros(3)) for r in history])
        return None

    def print_report(self):
        for r in self.results:
            print(f"\n{'='*60}")
            print(f"  Validation: {r['case_name']}")
            print(f"{'='*60}")
            print(f"  Simulation: T={r['T']}s, dt={r['dt']}s, steps={r['steps']}")
            for sname, sdata in r['signals'].items():
                m = sdata['metrics']
                print(f"  Signal: {sname}")
                print(f"    Ref final: {sdata['ref_last']:.6f}")
                print(f"    Our final: {sdata['our_last']:.6f}")
                print(f"    End error: {m['end_error']:.6f}")
                print(f"    RMSE:      {m['rmse']:.6f}")
                print(f"    Max error: {m['max_abs_error']:.6f}")
                print(f"    Status:    {sdata['status']}")
