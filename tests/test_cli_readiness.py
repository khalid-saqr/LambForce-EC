from pathlib import Path

from lambforce_ec.cli import main


ROOT = Path(__file__).resolve().parents[1]


def test_phase0_cli_performs_strict_repository_audit():
    assert main(["phase0-check", "--repository-root", str(ROOT)]) == 0
