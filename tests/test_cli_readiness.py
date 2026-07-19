from pathlib import Path

from lambforce_ec.cli import build_parser, main

ROOT = Path(__file__).resolve().parents[1]


def test_phase0_cli_performs_strict_repository_audit():
    assert main(["phase0-check", "--repository-root", str(ROOT)]) == 0


def test_cli_exposes_published_source_workflow_only():
    parser = build_parser()
    help_text = parser.format_help()
    choices = parser._subparsers._group_actions[0].choices
    assert "reproduce-hydrodynamics" in choices
    assert "verify-hydrodynamics" in choices
    assert "ingest-archive" not in choices
    assert "archive" not in help_text.lower()
