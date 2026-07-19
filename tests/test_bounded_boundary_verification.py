import numpy as np
from scipy import sparse

from lambforce_ec.solvers.bounded_fd import BoundedFDPlateSolver


def _biharmonic(n: int, length_x: float, length_z: float, support: str):
    dx = length_x / (n + 1)
    dz = length_z / (n + 1)
    d2x = BoundedFDPlateSolver._d2(n, dx)
    d2z = BoundedFDPlateSolver._d2(n, dz)
    d4x = BoundedFDPlateSolver._d4(n, dx, support)
    d4z = BoundedFDPlateSolver._d4(n, dz, support)
    ix = sparse.eye(n, format="csr")
    iz = sparse.eye(n, format="csr")
    return (
        sparse.kron(iz, d4x, format="csr")
        + 2 * sparse.kron(d2z, d2x, format="csr")
        + sparse.kron(d4z, ix, format="csr")
    )


def _manufactured_error(n: int, support: str) -> float:
    length_x = 1.0
    length_z = 0.8
    dx = length_x / (n + 1)
    dz = length_z / (n + 1)
    x = np.arange(1, n + 1) * dx
    z = np.arange(1, n + 1) * dz
    if support == "compliant_edge":
        ax = np.sin(np.pi * x / length_x)
        bz = np.sin(np.pi * z / length_z)
        field = ax[:, None] * bz[None, :]
        exact = ((np.pi / length_x) ** 2 + (np.pi / length_z) ** 2) ** 2 * field
    else:
        ax = x**2 * (length_x - x) ** 2
        bz = z**2 * (length_z - z) ** 2
        ax2 = 2 * length_x**2 - 12 * length_x * x + 12 * x**2
        bz2 = 2 * length_z**2 - 12 * length_z * z + 12 * z**2
        field = ax[:, None] * bz[None, :]
        exact = 24 * bz[None, :] + 2 * ax2[:, None] * bz2[None, :] + 24 * ax[:, None]
    approx = (
        _biharmonic(n, length_x, length_z, support)
        @ field.reshape(-1, order="F")
    ).reshape((n, n), order="F")
    return float(np.linalg.norm(approx - exact) / np.linalg.norm(exact))


def test_manufactured_support_solutions_converge():
    for support in ("compliant_edge", "clamped_edge"):
        coarse = _manufactured_error(15, support)
        fine = _manufactured_error(31, support)
        assert fine < coarse / 3.5
        assert fine < 2.0e-3


def test_clamped_fourth_derivative_is_exact_for_quartic_field():
    n = 31
    length = 1.0
    h = length / (n + 1)
    x = np.arange(1, n + 1) * h
    field = x**2 * (length - x) ** 2
    derivative = BoundedFDPlateSolver._d4(n, h, "clamped_edge") @ field
    assert np.max(np.abs(derivative - 24.0)) < 1e-9
