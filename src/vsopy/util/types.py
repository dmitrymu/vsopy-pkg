import astropy.units as u # type: ignore
from typing import NamedTuple

class ValErr(NamedTuple):
    """Value with uncertainty."""
    val: float
    """Value of the measurement."""
    err: float
    """Uncertainty of the measurement."""

ValErrDtype = [('val', 'f4'), ('err', 'f4')]

class MagErr(NamedTuple):
    """Magnitude with uncertainty."""
    mag: u.Quantity[u.mag]
    """Magnitude of the measurement."""
    err: u.Quantity[u.mag]
    """Uncertainty of the magnitude measurement."""

MagErrDtype = [('mag', 'f4'), ('err', 'f4')]
