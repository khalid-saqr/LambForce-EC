from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import numpy as np
import yaml

from .models import ArteryRecord, Geometry, MaterialState, RunConfig, SLSMaterial
from .structural import validate_correlated_glycocalyx
from .validation import checksum_arrays, require_keys, sha256_file
from .workflow import SimulationResult


def load_mapping(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream) if path.suffix.lower() in {".yaml", ".yml"} else json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a mapping.")
    return value


def load_artery_npz(path: str | Path) -> ArteryRecord:
    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        require_keys(
            set(data.files) and {name: None for name in data.files},
            {
                "radial_coordinate_m",
                "time_s",
                "lamb_density_signed_n_m3",
                "wall_shear_stress_pa",
                "metadata_json",
            },
            "artery NPZ",
        )
        metadata = json.loads(str(data["metadata_json"].item()))
        require_keys(
            metadata,
            {
                "artery_id",
                "artery_name",
                "radius_m",
                "omega0_rad_s",
                "source_identifier",
                "source_version",
                "source_checksum",
                "coordinate_convention",
            },
            "artery metadata",
        )
        iso = data["lamb_density_isotropic_n_m3"] if "lamb_density_isotropic_n_m3" in data else None
        record = ArteryRecord(
            artery_id=metadata["artery_id"],
            artery_name=metadata["artery_name"],
            radius_m=float(metadata["radius_m"]),
            omega0_rad_s=float(metadata["omega0_rad_s"]),
            radial_coordinate_m=data["radial_coordinate_m"],
            time_s=data["time_s"],
            lamb_density_signed_n_m3=data["lamb_density_signed_n_m3"],
            wall_shear_stress_pa=data["wall_shear_stress_pa"],
            lamb_density_isotropic_n_m3=iso,
            source_identifier=metadata["source_identifier"],
            source_version=metadata["source_version"],
            source_checksum=metadata["source_checksum"],
            coordinate_convention=metadata["coordinate_convention"],
            metadata=metadata.get("metadata", {}),
        )
    record.validate()
    return record


def save_artery_npz(record: ArteryRecord, path: str | Path) -> Path:
    record.validate()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "radial_coordinate_m": np.asarray(record.radial_coordinate_m),
        "time_s": np.asarray(record.time_s),
        "lamb_density_signed_n_m3": np.asarray(record.lamb_density_signed_n_m3),
        "wall_shear_stress_pa": np.asarray(record.wall_shear_stress_pa),
    }
    if record.lamb_density_isotropic_n_m3 is not None:
        payload["lamb_density_isotropic_n_m3"] = np.asarray(record.lamb_density_isotropic_n_m3)
    metadata = {
        "artery_id": record.artery_id,
        "artery_name": record.artery_name,
        "radius_m": record.radius_m,
        "omega0_rad_s": record.omega0_rad_s,
        "source_identifier": record.source_identifier,
        "source_version": record.source_version,
        "source_checksum": record.source_checksum,
        "coordinate_convention": record.coordinate_convention,
        "record_payload_checksum": checksum_arrays(payload, record.metadata),
        "metadata": record.metadata,
    }
    np.savez_compressed(path, **payload, metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)))
    return path


def config_from_mapping(mapping: dict[str, Any]) -> tuple[Geometry, MaterialState, RunConfig]:
    g = mapping.get("geometry", {})
    m = mapping.get("material", {})
    r = mapping.get("run", {})
    geometry = Geometry(**g)

    def sls(name: str, default: SLSMaterial) -> SLSMaterial:
        value = m.get(name)
        return default if value is None else SLSMaterial(**value)

    default_material = MaterialState()
    material = MaterialState(
        cortex=sls("cortex", default_material.cortex),
        cytosol=sls("cytosol", default_material.cytosol),
        glycocalyx=sls("glycocalyx", default_material.glycocalyx),
        nucleus=sls("nucleus", default_material.nucleus),
        poisson_ratio=float(m.get("poisson_ratio", default_material.poisson_ratio)),
        parameter_set_id=str(m.get("parameter_set_id", default_material.parameter_set_id)),
        glycocalyx_state_id=str(m.get("glycocalyx_state_id", default_material.glycocalyx_state_id)),
    )
    config = RunConfig(**r)
    geometry.validate()
    material.validate()
    config.validate()
    validate_correlated_glycocalyx(geometry, material)
    return geometry, material, config


def save_result(result: SimulationResult, output_directory: str | Path) -> Path:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    npz_path = output / "fields.npz"
    np.savez_compressed(npz_path, **result.arrays)
    units = {
        "load_normal_pa": "Pa",
        "wall_shear_stress_pa": "Pa",
        "displacement_normal_cell_m": "m",
        "displacement_normal_apical_top_m": "m",
        "displacement_tangential_m": "m",
        "curvature_x_m_inv": "m^-1",
        "curvature_z_m_inv": "m^-1",
        "curvature_xz_m_inv": "m^-1",
        "curvature_change_m_inv": "m^-1",
        "tension_x_n_m": "N m^-1",
        "tension_z_n_m": "N m^-1",
        "tension_xz_total_n_m": "N m^-1",
        "membrane_tension_max_principal_n_m": "N m^-1",
        "strain_max_principal": "1",
        "glycocalyx_strain_normal": "1",
        "glycocalyx_reaction_stress_pa": "Pa",
        "strain_energy_density_j_m3": "J m^-3",
        "tension_loading_rate_n_m_s": "N m^-1 s^-1",
    }
    manifest = {
        "metadata": result.metadata,
        "summaries": result.summaries,
        "fields_file": npz_path.name,
        "fields_sha256": sha256_file(npz_path),
        "array_units": {name: unit for name, unit in units.items() if name in result.arrays},
    }
    with (output / "manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
    return output
