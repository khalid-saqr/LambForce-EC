from __future__ import annotations

import argparse
import json
from pathlib import Path

from .archive import ingest_archive_member, qualify_hydrodynamics, save_qualification_report
from .io import (
    config_from_mapping,
    load_artery_npz,
    load_mapping,
    save_artery_npz,
    save_result,
)
from .protocol import (
    load_yaml,
    resolve_data_path,
    validate_parameter_registry,
    validate_reference_arteries,
    validate_source_registry,
    validate_traceability_matrix,
)
from .synthetic import make_synthetic_artery
from .validation import sha256_file
from .workflow import run_case


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lambforce-ec")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="Validate an artery NPZ and its payload checksum.")
    validate.add_argument("input")

    synthetic = sub.add_parser("synthetic", help="Create a non-claim-bearing software-validation input.")
    synthetic.add_argument("output")

    run = sub.add_parser("run", help="Run one registered mechanics case.")
    run.add_argument("--input", required=True)
    run.add_argument("--config", required=True)
    run.add_argument("--output", required=True)
    run.add_argument(
        "--source-registry",
        help="Versioned source registry; required for claim-bearing execution.",
    )

    phase0 = sub.add_parser("phase0-check", help="Validate readiness registries and traceability.")
    phase0.add_argument("--parameters")
    phase0.add_argument("--traceability")
    phase0.add_argument("--source-registry")
    phase0.add_argument("--reference-arteries")
    phase0.add_argument(
        "--repository-root",
        help="Repository root for strict path/test/README traceability checks.",
    )

    ingest = sub.add_parser(
        "ingest-archive",
        help="Deterministically convert one immutable hydrodynamic archive member.",
    )
    ingest.add_argument("--archive", required=True)
    ingest.add_argument("--member-npz", required=True)
    ingest.add_argument("--manifest", required=True)
    ingest.add_argument("--output", required=True)
    ingest.add_argument("--converter-commit", required=True)

    qualify = sub.add_parser(
        "verify-hydrodynamics",
        help="Verify one converted record against the registered immutable archive and harmonics.",
    )
    qualify.add_argument("--input", required=True)
    qualify.add_argument("--source-registry")
    qualify.add_argument("--output", required=True)
    return parser


def _repository_root(value: str | None) -> Path | None:
    if value is not None:
        root = Path(value).resolve()
        if not (root / "README.md").is_file():
            raise ValueError(f"{root} is not a repository root containing README.md.")
        return root
    current = Path.cwd().resolve()
    if (current / "README.md").is_file() and (current / "src/lambforce_ec").is_dir():
        return current
    return None


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "validate":
        record = load_artery_npz(args.input)
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "artery_id": record.artery_id,
                    "record_payload_sha256": record.record_payload_sha256,
                },
                indent=2,
            )
        )
        return 0
    if args.command == "synthetic":
        path = save_artery_npz(make_synthetic_artery(), args.output)
        print(path)
        return 0
    if args.command == "run":
        record = load_artery_npz(args.input)
        geometry, material, config = config_from_mapping(load_mapping(args.config))
        source_registry = None
        registry_path = None
        if config.claim_bearing:
            registry_path = resolve_data_path(
                "registry/source_registry.yaml", args.source_registry
            )
            source_registry = load_yaml(registry_path)
        result = run_case(
            record,
            geometry,
            material,
            config,
            source_registry=source_registry,
        )
        result.metadata.update(
            {
                "record_payload_sha256": record.record_payload_sha256,
                "source_member_sha256": record.source_member_sha256,
                "conversion_manifest_sha256": record.conversion_manifest_sha256,
                "converter_commit_sha": record.converter_commit_sha,
            }
        )
        if registry_path is not None:
            result.metadata["source_registry_sha256"] = sha256_file(registry_path)
        output = save_result(result, args.output)
        print(output)
        return 0
    if args.command == "phase0-check":
        parameter_path = resolve_data_path(
            "registry/parameter_registry.csv", args.parameters
        )
        traceability_path = resolve_data_path(
            "protocol/readme_traceability.yaml", args.traceability
        )
        source_path = resolve_data_path(
            "registry/source_registry.yaml", args.source_registry
        )
        reference_path = resolve_data_path(
            "configs/reference_arteries.yaml", args.reference_arteries
        )
        rows = validate_parameter_registry(parameter_path)
        source_count = validate_source_registry(load_yaml(source_path))
        artery_count = validate_reference_arteries(load_yaml(reference_path))
        counts = validate_traceability_matrix(
            load_yaml(traceability_path),
            repository_root=_repository_root(args.repository_root),
        )
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "parameter_rows": rows,
                    "source_records": source_count,
                    "reference_arteries": artery_count,
                    **counts,
                },
                indent=2,
            )
        )
        return 0
    if args.command == "ingest-archive":
        report = ingest_archive_member(
            args.archive,
            args.member_npz,
            args.manifest,
            args.output,
            args.converter_commit,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "verify-hydrodynamics":
        record = load_artery_npz(args.input)
        registry_path = resolve_data_path(
            "registry/source_registry.yaml", args.source_registry
        )
        report = qualify_hydrodynamics(record, load_yaml(registry_path))
        save_qualification_report(report, args.output)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    raise RuntimeError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
