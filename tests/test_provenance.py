import pytest

from lambforce_ec.exceptions import ProvenanceError
from lambforce_ec.models import RunConfig
from lambforce_ec.provenance import compute_record_payload_sha256
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case


def test_synthetic_record_cannot_be_promoted_to_claim_bearing():
    record = make_synthetic_artery(n_time=16)
    record.record_payload_sha256 = compute_record_payload_sha256(record)
    config = RunConfig(nx=8, nz=8, claim_bearing=True)
    forged_registry = {
        "registry_version": "forged",
        "sources": [
            {
                "source_identifier": record.source_identifier,
                "source_version": record.source_version,
                "artery_ids": [record.artery_id],
                "source_status": "synthetic_non_claim_bearing",
                "repository_commit_sha": None,
                "published_notebook_path": None,
                "published_notebook_blob_sha": None,
                "published_input_registry_sha256": None,
                "reproduction_commit_sha": None,
                "record_payload_sha256_by_mode_and_artery": {},
                "qualification_status": "not_applicable",
                "role": "software_validation_only",
                "claim_bearing": False,
            }
        ],
    }
    with pytest.raises(ProvenanceError):
        run_case(record, config=config, source_registry=forged_registry)
