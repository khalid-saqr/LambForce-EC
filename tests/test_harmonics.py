import numpy as np

from lambforce_ec.harmonics import apply_harmonic_control, forward_real_series, inverse_real_series


def test_fft_round_trip_and_registered_controls():
    n = 128
    period = 0.8
    t = np.arange(n) * period / n
    x = 1.2 + np.sin(2 * np.pi * t / period) + 0.4 * np.sin(4 * np.pi * t / period + 0.3)
    series = forward_real_series(x, t)
    assert np.max(np.abs(inverse_real_series(series) - x)) < 1e-12
    fundamental = apply_harmonic_control(series, "fundamental_only")
    assert np.all(fundamental.coefficients[2:] == 0)
    le2 = apply_harmonic_control(series, "harmonics_le_2")
    assert np.all(le2.coefficients[3:] == 0)
