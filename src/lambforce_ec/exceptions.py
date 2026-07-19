class LambForceECError(Exception):
    """Base package error."""


class ValidationError(LambForceECError):
    """Raised when a frozen contract or input validation rule is violated."""


class ConservationError(LambForceECError):
    """Raised when applied and reaction resultants do not agree."""


class ProvenanceError(LambForceECError):
    """Raised when a claim-bearing source cannot be resolved to a verified archive."""
