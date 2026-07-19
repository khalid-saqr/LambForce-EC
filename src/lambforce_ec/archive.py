"""Deprecated archive-import compatibility surface.

The publication source is an exact Git notebook artifact, not a presumed external artery-member
archive. Use :mod:`lambforce_ec.published_source` and the `reproduce-hydrodynamics` CLI command.
"""

from __future__ import annotations

from .exceptions import ValidationError


def ingest_archive_member(*args, **kwargs):
    del args, kwargs
    raise ValidationError(
        "ingest-archive was removed in version 0.6.0. Use reproduce-hydrodynamics with the "
        "frozen published_v2_hydrodynamics.yaml registry."
    )


def qualify_hydrodynamics(*args, **kwargs):
    del args, kwargs
    raise ValidationError(
        "Single archive-member qualification is obsolete. Use verify-hydrodynamics on the "
        "complete all-six published-source reproduction directory."
    )
