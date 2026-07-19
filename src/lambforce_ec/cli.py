from __future__ import annotations

import argparse
import json

from .io import config_from_mapping, load_artery_npz, load_mapping, save_artery_npz, save_result
from .protocol import load_yaml, validate_parameter_registry, validate_traceability_matrix
from .synthetic import make_synthetic_artery
from .workflow import run_case


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lambforce-ec")
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate", help="Validate an artery NPZ against the input contract.")
    validate.add_argument("input")
    synthetic = sub.add_parser("synthetic", help="Create a non-claim-bearing software-validation input.")
    synthetic.add_argument("output")
    run = sub.add_parser("run", help="Run one registered mechanics case.")
    run.add_argument("--input", required=True)
    run.add_argument("--config", required=True)
    run.add_argument("--output", required=True)
    phase0 = sub.add_parser("phase0-check", help="Validate Phase 0 registries and traceability.")
    phase0.add_argument("--parameters", default="registry/parameter_registry.csv")
    phase0.add_argument("--traceability", default="protocol/readme_traceability.yaml")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "validate":
        record = load_artery_npz(args.input)
        print(json.dumps({"status": "PASS", "artery_id": record.artery_id}, indent=2))
        return 0
    if args.command == "synthetic":
        path = save_artery_npz(make_synthetic_artery(), args.output)
        print(path)
        return 0
    if args.command == "run":
        record = load_artery_npz(args.input)
        geometry, material, config = config_from_mapping(load_mapping(args.config))
        output = save_result(run_case(record, geometry, material, config), args.output)
        print(output)
        return 0
    if args.command == "phase0-check":
        rows = validate_parameter_registry(args.parameters)
        counts = validate_traceability_matrix(load_yaml(args.traceability))
        print(json.dumps({"status": "PASS", "parameter_rows": rows, **counts}, indent=2))
        return 0
    raise RuntimeError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
