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
                "archive_status": "verified",
                "archive_sha256": record.source_archive_sha256,
                "source_member_sha256_by_artery": {
                    record.artery_id: record.source_member_sha256
                },
                "conversion_manifest_sha256_by_artery": {
                    record.artery_id: record.conversion_manifest_sha256
                },
                "converter_commit_sha": record.converter_commit_sha,
                "role": "software_validation_only",
                "claim_bearing": False,
            }
        ],
    }
    with pytest.raises(ProvenanceError):
        run_case(record, config=config, source_registry=forged_registry)
