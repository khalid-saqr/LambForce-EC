class LambForceECError(Exception):
    """Base package exception."""


class ValidationError(LambForceECError, ValueError):
    """Raised when a scientific input record violates the data contract."""


class ConservationError(LambForceECError, RuntimeError):
    """Raised when a load or reaction resultant fails the registered tolerance."""
