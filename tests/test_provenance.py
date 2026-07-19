import copy
import pytest

from lambforce_ec.exceptions import ProvenanceError
from lambforce_ec.models import RunConfig
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case


def test_claim_bearing_requires_verified_archive():
    record = make_synthetic_artery(n_time=16)
    config = RunConfig(nx=8, nz=8, claim_bearing=True)
    with pytest.raises(ProvenanceError):
        run_case(record, config=config)
    registry = {
        "sources": [{
            "source_identifier": record.source_identifier,
            "source_version": record.source_version,
            "artery_ids": [record.artery_id],
            "archive_status": "awaiting_archive",
            "archive_sha256": record.source_checksum,
        }]
    }
    with pytest.raises(ProvenanceError):
        run_case(record, config=config, source_registry=registry)
    verified = copy.deepcopy(registry)
    verified["sources"][0]["archive_status"] = "verified"
    result = run_case(record, config=config, source_registry=verified)
    assert result.metadata["claim_bearing"] is True
