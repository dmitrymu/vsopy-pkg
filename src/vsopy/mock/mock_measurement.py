import astropy.units as u
import numpy as np
from astropy.table import QTable

def mock_measure(stars:QTable, band:tuple[str, str],
                 Rab:tuple[float, float, float], Ra:tuple[float, float, float]) -> QTable:
    bandA, bandB = band
    N = len(stars)
    A = stars[bandA]['mag']
    B = stars[bandB]['mag']
    C = A - B
    Tab, Zab, Sab = Rab
    c = C * (1/Tab) + Zab * u.mag + np.random.normal(0, Sab, size=N) * u.mag

    Ta, Za, Sa = Ra
    a = A - C * Ta + Za * u.mag + np.random.normal(0, Sa, size=N) * u.mag
    b = a - c

    return QTable({
        'auid': stars['auid'],
        bandA: a,
        bandB: b
    })
