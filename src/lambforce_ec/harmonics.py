from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .exceptions import ValidationError


@dataclass(frozen=True)
class HarmonicSeries:
    coefficients: np.ndarray
    omega_rad_s: np.ndarray
    n_time: int
    dt_s: float


def forward_real_series(values: np.ndarray, time_s: np.ndarray) -> HarmonicSeries:
    x = np.asarray(values, dtype=float)
    t = np.asarray(time_s, dtype=float)
    if x.shape[0] != t.size or t.ndim != 1 or t.size < 4:
        raise ValidationError("time series leading dimension must match time_s.")
    dt = np.diff(t)
    if not np.allclose(dt, dt[0], rtol=1e-8, atol=1e-12):
        raise ValidationError("time_s must be uniformly sampled.")
    coefficients = np.fft.rfft(x, axis=0) / t.size
    frequencies_hz = np.fft.rfftfreq(t.size, d=float(dt[0]))
    return HarmonicSeries(coefficients, 2 * np.pi * frequencies_hz, t.size, float(dt[0]))


def inverse_real_series(series: HarmonicSeries) -> np.ndarray:
    return np.fft.irfft(series.coefficients * series.n_time, n=series.n_time, axis=0)


def apply_harmonic_control(series: HarmonicSeries, control: str) -> HarmonicSeries:
    coeff = np.array(series.coefficients, copy=True)
    if control == "full_waveform":
        pass
    elif control == "fundamental_only":
        if coeff.shape[0] > 2:
            coeff[2:] = 0
    elif control == "harmonics_le_2":
        if coeff.shape[0] > 3:
            coeff[3:] = 0
    else:
        raise ValidationError(f"Unknown harmonic control: {control}")
    return HarmonicSeries(coeff, series.omega_rad_s, series.n_time, series.dt_s)


def reconstruct_coefficients(coefficients: np.ndarray, n_time: int) -> np.ndarray:
    return np.fft.irfft(np.asarray(coefficients) * n_time, n=n_time, axis=0)
