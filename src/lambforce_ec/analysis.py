from __future__ import annotations

from typing import Any
import numpy as np


def principal_max(a: np.ndarray, b: np.ndarray, shear: np.ndarray) -> np.ndarray:
    center = 0.5 * (a + b)
    radius = np.sqrt((0.5 * (a - b)) ** 2 + shear**2)
    return center + radius


def field_summary(field: np.ndarray, time_s: np.ndarray) -> dict[str, Any]:
    x = np.asarray(field, dtype=float)
    t = np.asarray(time_s, dtype=float)
    abs_x = np.abs(x)
    flat_index = int(np.argmax(abs_x))
    index = np.unravel_index(flat_index, x.shape)
    time_index = index[0] if x.ndim >= 1 and x.shape[0] == t.size else None
    spatial_abs_mean = np.mean(abs_x, axis=tuple(range(1, x.ndim))) if x.ndim > 1 else abs_x
    peak_to_mean = float(np.max(abs_x) / max(float(np.mean(abs_x)), 1e-30))
    return {
        "peak_signed": float(x[index]),
        "peak_absolute": float(np.max(abs_x)),
        "rms": float(np.sqrt(np.mean(x**2))),
        "cycle_mean": float(np.mean(x)),
        "percentile_95_absolute": float(np.percentile(abs_x, 95)),
        "peak_to_peak": float(np.max(x) - np.min(x)),
        "time_of_peak_s": None if time_index is None else float(t[time_index]),
        "index_of_peak": [int(i) for i in index],
        "spatial_peak_to_mean": peak_to_mean,
        "cycle_mean_spatial_absolute": float(np.mean(spatial_abs_mean)),
    }


def incremental_field(combined: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    a = np.asarray(combined)
    b = np.asarray(baseline)
    if a.shape != b.shape:
        raise ValueError("combined and baseline fields must have the same shape.")
    return a - b


def norm_ratio(numerator: np.ndarray, denominator: np.ndarray, floor: float = 1e-30) -> float:
    n = float(np.linalg.norm(np.asarray(numerator).ravel()))
    d = float(np.linalg.norm(np.asarray(denominator).ravel()))
    return float(n / max(d, floor))
