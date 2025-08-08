import astropy.units as u
import numpy as np
import unittest

from astropy.coordinates import SkyCoord
from astropy.table import QTable
from vsopy.mock import MockImageBuilder, MockStar
from vsopy.phot import measure_photometry, filter_centroids
from vsopy.util import Aperture

SHAPE = (31,31)
NOISE_MEAN = 100
NOISE_STDDEV = 10
STAR_PEAK = 1000
STAR_POS = (15, 15)
STAR_FWHM = 3
STAR_AUID = 'mock-star-1'
STAR_AUID_2 = 'mock-star-2'
STAR_AUID_3 = 'mock-star-3'

class MeasurePhotometryTest(unittest.TestCase):

    def test_one_star(self):

        builder = MockImageBuilder(SHAPE)
        builder.add_noise(NOISE_MEAN, NOISE_STDDEV)
        builder.add_star(MockStar(STAR_PEAK, STAR_POS, STAR_FWHM, 0, 0*u.deg))
        image = builder.get_image(1.2 * u.arcsec,
                        dict(
                            EXPTIME=2,
                            GAIN=100,
                            instrume="ZWO CCD ASI533MM Pro"
                        ))
        centroids = QTable(dict(
            auid = [STAR_AUID, STAR_AUID_2],
            radec2000 = SkyCoord(ra=[0, 3600] * u.arcsec, dec=[-0, 3600] * u.arcsec)
        ))

        ph = measure_photometry(image, centroids, Aperture(5, 10, 15))

        self.assertSequenceEqual(
            ph.colnames, ['auid', 'radec2000', 'flux', 'snr', 'M', 'peak', 'sky_centroid'])
        self.assertEqual(len(ph), 1)
        self.assertEqual(ph['auid'][0], STAR_AUID)

    def test_filter_centroids(self):

        builder = MockImageBuilder(SHAPE)
        builder.add_noise(NOISE_MEAN, NOISE_STDDEV)
        builder.add_star(MockStar(STAR_PEAK, STAR_POS, STAR_FWHM, 0, 0*u.deg))
        image = builder.get_image(1.2 * u.arcsec,
                        dict(
                            EXPTIME=2,
                            GAIN=100,
                            instrume="ZWO CCD ASI533MM Pro"
                        ))
        # stars: inside wih aperture, outside with aperture, inside with partial aperture
        centroids = QTable(dict(
            auid = [STAR_AUID, STAR_AUID_2, STAR_AUID_3],
            radec2000 = SkyCoord(ra=[0, 3600, 16] * u.arcsec, dec=[-0, 3600, 16] * u.arcsec)
        ))

        filtered = filter_centroids(image, centroids, 5 * u.arcsec)

        self.assertSequenceEqual(filtered.colnames, ['auid', 'radec2000'])
        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered['auid'][0], STAR_AUID)
