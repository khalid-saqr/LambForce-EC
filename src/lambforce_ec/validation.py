from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np

from .exceptions import ConservationError, ValidationError


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def checksum_arrays(arrays: dict[str, np.ndarray], metadata: dict[str, Any] | None = None) -> str:
    h = hashlib.sha256()
    for key in sorted(arrays):
        array = np.ascontiguousarray(arrays[key])
        h.update(key.encode("utf-8"))
        h.update(str(array.dtype).encode("ascii"))
        h.update(canonical_json_bytes(list(array.shape)))
        h.update(array.tobytes())
    if metadata is not None:
        h.update(canonical_json_bytes(metadata))
    return h.hexdigest()


def relative_error(value: complex | float, reference: complex | float, floor: float = 1e-30) -> float:
    return float(abs(value - reference) / max(abs(reference), floor))


def assert_conserved(
    applied: complex | float,
    reaction: complex | float,
    tolerance_relative: float,
    label: str,
) -> None:
    error = relative_error(reaction, applied)
    if error > tolerance_relative:
        raise ConservationError(
            f"{label} resultant conservation failed: relative error={error:.3e}, "
            f"tolerance={tolerance_relative:.3e}."
        )


def require_keys(mapping: dict[str, Any], required: set[str], context: str) -> None:
    missing = sorted(required - set(mapping))
    if missing:
        raise ValidationError(f"{context} is missing required keys: {', '.join(missing)}")
