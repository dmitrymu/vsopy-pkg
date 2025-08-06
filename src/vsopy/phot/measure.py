import astropy.units as u # type: ignore
import numpy as np
from .. import reduce
from ..data import CameraRegistry
from astropy.nddata import CCDData # type: ignore
from astropy.table import QTable, Column # type: ignore
from photutils.aperture import (ApertureStats, # type: ignore
                                SkyCircularAperture,
                                SkyCircularAnnulus)
from vsopy.util import Aperture


def measure_photometry(image:CCDData, stars:QTable, aperture:Aperture,
                       extended:bool=False) -> QTable:
    """Extract aperture photometry data from the image.

    Given the calibrated image and the list of star centroids in sky coordinates,
    measure star flux using circular aperture and sky annulus.  For detailed discussion
    see Warner B.D., A Practical Guide to Lightcurve Photometry and Analysis, 4.2.

    Flux counting is performed by photutils package
    (https://photutils.readthedocs.io). Uncertainty fromthe image is propagated to
    the flux values. Flux is normalized by the image exposure to electrons per second.

    The flux for central aperture :math:`F_c` contains both star and sky electrons, the flux
    for annulus :math:`F_a` - sky electrons only. Number of sky electrons in the central
    aperture is

    .. math::
        :label: measure.1

        F_{sky} = \\frac{F_{a}}{N_{ann}} N_c,

    where :math:`N_{ann}` and :math:`N_c` are pixel counts for the annulus
    and the central aperture respectively.

    Then star flux in electrons is

    .. math::
        :label: apply.2

        F_{star} = F_c - F_{sky}

    and signal to noise ratio is

    .. math::
        :label: apply.3

        SNR = \\frac{F_{star}}{\\sqrt{F_{star} + 2 F_{sky}}}

    To convert flux to magnitude, an arbitrary scale value
    :math:`F_0 = 10^6 \\frac{e}{s}` is used:

    .. math::
        :label: apply.4

        M = -2.5 log_{10}(F_{star} / F_0)

    :param image:       calibrated image, pixel counts in electrons.
    :type image:        :py:class:`~astropy.nddata.CCDData`
    :param stars:       list of stars to be measured: AUID, centroid, optional name
    :type stars:        :py:class:`~astropy.table.QTable`
    :param aperture:    circular aperture and annulus radii in sky coordinates
    :type aperture:     Aperture
    :param extended:    whether to return additional stats (FWHM,
                        ellipticity etc), default False
    :type extended: bool
    :return:           photometry results
    :rtype:            :py:class:`~astropy.table.QTable`, fields:

                - auid,
                - centroid,
                - flux (:math:`\\frac{e}{s}`)
                - SNR (db)
                - Magnitude with uncertainty,
                - relative peak (max count / saturation level)
                - ellipticity (if extended=True)
                - FWHM (if extended=True)
                - orientation (if extended=True)
    """
    result = QTable(stars['auid', 'radec2000'])
    if 'name' in stars.colnames:
        result['name'] = stars['name']

    corrected = CCDData(image)

    apr = SkyCircularAperture(stars['radec2000'], r=aperture.r)
    ap_stats = ApertureStats(corrected, apr)

    ann = SkyCircularAnnulus(stars['radec2000'],
                             r_in=aperture.r_in,
                             r_out=aperture.r_out)
    ann_stats = ApertureStats(corrected, ann)

    zero_level = 1_000_000 * image.unit / u.second
    exp = image.header['EXPTIME'] * u.second
    camera = CameraRegistry.get(image.header['instrume'])
    read_noise = 0 if not camera or image.unit != u.electron else camera.read_noise(image.header['GAIN'])
    Ft = ap_stats.sum
    Fb = ap_stats.sum_aper_area.value * ann_stats.mean
    FtFb = Ft - Fb
    Ft2Fb = Ft + 2 * Fb + read_noise * ap_stats.sum_aper_area.value
    result['flux'] = FtFb / exp
    # prevent NaN in snr column
    snr = np.clip(FtFb.value/np.sqrt(np.abs(Ft2Fb.value)), 1e-10, 1e30)
    result['snr'] = 10 * np.log10(snr) * u.db
    mag = -2.5 * np.log10(result['flux'] / zero_level).value

    Ft_err = ap_stats.sum_err
    Fb_err = ap_stats.sum_aper_area.value * (ann_stats.sum_err / ann_stats.sum_aper_area.value)
    flux_err = np.sqrt(Ft_err*Ft_err + Fb_err*Fb_err) / exp
    mag_err = (2.5 * flux_err / result['flux'] / np.log(10)).value
    result['M'] = Column(list(zip(mag, mag_err)),
                         unit=u.mag,
                         dtype=[('mag','f4'), ('err', 'f4')])

    result['peak'] = np.nan if not camera else ap_stats.max.value / camera.max_adu.value
    result['sky_centroid'] = ap_stats.sky_centroid

    if extended:
        result['ellipticity'] = ap_stats.ellipticity
        result['fwhm'] = ap_stats.fwhm
        result['orientation'] = ap_stats.orientation

    return result


def process_image(path, matcher, solver, centroids, aperture):
    try:
        image = reduce.update_wcs(CCDData.read(path, unit='adu'), solver(path))
        calibration = matcher.match(image.header)
        reduced = reduce.calibrate_image(image,
                                        dark=calibration.dark,
                                        flat=calibration.flat)
        return measure_photometry(reduced, centroids(reduced), aperture)
    except Exception:
        return None
