from __future__ import annotations

import argparse
import json
from pathlib import Path

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
from .published_source import (
    reproduce_all_six,
    save_reproduction_verification,
    validate_published_inputs,
    validate_published_source_binding,
    verify_reproduction_directory,
)
from .synthetic import make_synthetic_artery
from .validation import sha256_file
from .workflow import run_case


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lambforce-ec")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="Validate an artery NPZ and its payload checksum.")
    validate.add_argument("input")

    synthetic = sub.add_parser(
        "synthetic", help="Create a non-claim-bearing software-validation input."
    )
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
    phase0.add_argument("--published-inputs")
    phase0.add_argument(
        "--repository-root",
        help="Repository root for strict path/test/README traceability checks.",
    )

    reproduce = sub.add_parser(
        "reproduce-hydrodynamics",
        help="Reproduce all six published hydrodynamic records from the frozen v2 source inputs.",
    )
    reproduce.add_argument("--published-inputs")
    reproduce.add_argument("--output", required=True)
    reproduce.add_argument("--reproduction-commit", required=True)
    reproduce.add_argument(
        "--profile",
        choices=("verification", "publication"),
        default="publication",
    )

    qualify = sub.add_parser(
        "verify-hydrodynamics",
        help="Verify an all-six published-source reproduction directory.",
    )
    qualify.add_argument("--reproduction-directory", required=True)
    qualify.add_argument("--published-inputs")
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
                    "published_source": record.metadata.get("published_source"),
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
            registry_path = resolve_data_path("registry/source_registry.yaml", args.source_registry)
            source_registry = load_yaml(registry_path)
        result = run_case(
            record,
            geometry,
            material,
            config,
            source_registry=source_registry,
        )
        published = record.metadata.get("published_source", {})
        result.metadata.update(
            {
                "record_payload_sha256": record.record_payload_sha256,
                "source_repository_commit_sha": published.get("repository_commit_sha"),
                "source_notebook_blob_sha": published.get("published_notebook_blob_sha"),
                "published_input_registry_sha256": published.get(
                    "published_input_registry_sha256"
                ),
                "reproduction_commit_sha": published.get("reproduction_commit_sha"),
                "reproduction_mode": published.get("reproduction_mode"),
            }
        )
        if registry_path is not None:
            result.metadata["source_registry_sha256"] = sha256_file(registry_path)
        output = save_result(result, args.output)
        print(output)
        return 0
    if args.command == "phase0-check":
        parameter_path = resolve_data_path("registry/parameter_registry.csv", args.parameters)
        traceability_path = resolve_data_path(
            "protocol/readme_traceability.yaml", args.traceability
        )
        source_path = resolve_data_path("registry/source_registry.yaml", args.source_registry)
        reference_path = resolve_data_path(
            "configs/reference_arteries.yaml", args.reference_arteries
        )
        published_path = resolve_data_path(
            "registry/published_v2_hydrodynamics.yaml", args.published_inputs
        )
        rows = validate_parameter_registry(parameter_path)
        source_registry = load_yaml(source_path)
        source_count = validate_source_registry(source_registry)
        artery_count = validate_reference_arteries(load_yaml(reference_path))
        published_inputs = load_yaml(published_path)
        published_count = validate_published_inputs(published_inputs)
        published_sha = sha256_file(published_path)
        validate_published_source_binding(published_inputs, published_sha, source_registry)
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
                    "published_source_arteries": published_count,
                    "published_input_registry_sha256": published_sha,
                    **counts,
                },
                indent=2,
            )
        )
        return 0
    if args.command == "reproduce-hydrodynamics":
        published_path = resolve_data_path(
            "registry/published_v2_hydrodynamics.yaml", args.published_inputs
        )
        report = reproduce_all_six(
            published_path,
            args.output,
            args.reproduction_commit,
            args.profile,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if args.command == "verify-hydrodynamics":
        published_path = resolve_data_path(
            "registry/published_v2_hydrodynamics.yaml", args.published_inputs
        )
        registry_path = resolve_data_path(
            "registry/source_registry.yaml", args.source_registry
        )
        report = verify_reproduction_directory(
            args.reproduction_directory,
            published_path,
            load_yaml(registry_path),
        )
        save_reproduction_verification(report, args.output)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    raise RuntimeError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
