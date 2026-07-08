import numpy as np


def compute_error_metrics(ref_signal, our_signal, ref_time, our_time):
    our_interp = np.interp(ref_time, our_time, our_signal)
    diff = our_interp - ref_signal
    rmse = np.sqrt(np.mean(diff ** 2))
    max_abs_error = np.max(np.abs(diff))
    mean_abs_error = np.mean(np.abs(diff))
    end_error = diff[-1] if len(diff) > 0 else 0.0

    ref_range = np.max(ref_signal) - np.min(ref_signal) if len(ref_signal) > 0 else 1.0
    normalized_rmse = rmse / ref_range if ref_range > 1e-30 else 0.0

    return {
        'rmse': float(rmse),
        'max_abs_error': float(max_abs_error),
        'mean_abs_error': float(mean_abs_error),
        'end_error': float(end_error),
        'normalized_rmse': float(normalized_rmse),
        'num_points': len(diff),
    }


def classify_error(metrics, smoke_tol=0.25, strict_tol=0.1):
    end_err = abs(metrics['end_error'])
    if end_err < strict_tol:
        return 'PASS (strict)', 'strict'
    elif end_err < smoke_tol:
        return 'PASS (smoke)', 'smoke'
    else:
        return 'FAIL', 'fail'
