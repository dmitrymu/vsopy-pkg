import astropy.units as u
import numpy as np
import unittest

from unittest.mock import patch, mock_open
from vsopy.util import Aperture, Settings

class ApertureTest(unittest.TestCase):

    def test_constructor(self):

        ap = Aperture(1, 2, 3)
        self.assertEqual(ap.r, 1*u.arcsec)
        self.assertEqual(ap.r_in, 2*u.arcsec)
        self.assertEqual(ap.r_out, 3*u.arcsec)

        ap = Aperture(1, 2, 3*u.arcmin)
        self.assertEqual(ap.r, 1*u.arcsec)
        self.assertEqual(ap.r_in, 2*u.arcsec)
        self.assertEqual(ap.r_out, 180*u.arcsec)

        ap = Aperture(1, 2, 3, unit=u.arcmin)
        self.assertEqual(ap.r, 1*u.arcmin)
        self.assertEqual(ap.r_in, 2*u.arcmin)
        self.assertEqual(ap.r_out, 3*u.arcmin)

    def test_conversion(self):

        ap= Aperture(1, 2, 3)
        d = ap.to_dict()
        self.assertDictEqual(d, dict(r_ap=1, r_in=2, r_out=3, unit=str(u.arcsec)))

        ap1 = Aperture.from_dict(d)
        self.assertEqual(ap1.r, 1*u.arcsec)
        self.assertEqual(ap1.r_in, 2*u.arcsec)
        self.assertEqual(ap1.r_out, 3*u.arcsec)

        ap2 = ap1.to_pixels(2*u.arcsec/u.pixel)
        self.assertEqual(ap2.r, .5*u.pixel)
        self.assertEqual(ap2.r_in, 1.0*u.pixel)
        self.assertEqual(ap2.r_out, 1.5*u.pixel)

class SettingsTest(unittest.TestCase):

    DEFAULT_JSON="""
{
  "aperture": {
    "unit": "arcsec",
    "r_ap": 2.0,
    "r_in": 4.0,
    "r_out": 8.0
  },
  "disabled": [
    "mock-star-1"
  ],
  "bands": [
    "X",
    "Y"
  ],
  "diff_photometry": {
    "XY": {
      "comp": "comp-star-1",
      "check": "check-star-2",
      "start": 0,
      "finish": 10,
      "rejected": [1, 2, 3]
    }
  }
}
"""

    def test_constructor(self):
        s = Settings(None)
        with self.assertRaises(KeyError):
            s.aperture
        self.assertSequenceEqual(s.bands, [])

    def test_disabled_star(self):
        s = Settings(None)
        self.assertTrue(s.is_star_enabled('star1'))
        self.assertTrue(s.is_star_enabled('star2'))
        self.assertTrue(s.is_star_enabled('star1'))
        self.assertTrue(s.is_star_enabled('star2'))

        s.disable_star('star1')
        self.assertFalse(s.is_star_enabled('star1'))
        self.assertTrue(s.is_star_enabled('star2'))

        s.disable_star('star1')  # assert no exception

        s.enable_star('star1')
        self.assertTrue(s.is_star_enabled('star1'))
        self.assertTrue(s.is_star_enabled('star2'))

        s.enable_star('star1') # assert no exception

    def test_photometry_settings(self):
        s = Settings(None)
        self.assertSetEqual(s.photometry(('X', 'Y')).rejected, set())

        s.photometry(('X', 'Y')).set_check('check-star-1')
        self.assertEqual(s.data_["diff_photometry"]["XY"]["check"], 'check-star-1')

        s.photometry(('X', 'Y')).set_rejected(None)
        self.assertSetEqual(s.photometry(('X', 'Y')).rejected, set())

        s.photometry(('X', 'Y')).set_rejected([])
        self.assertSetEqual(s.photometry(('X', 'Y')).rejected, set())

        s.photometry(('X', 'Y')).set_rejected(np.array([42]))
        self.assertSetEqual(s.photometry(('X', 'Y')).rejected, set([42]))

    @patch('vsopy.util.Settings.Path.exists', return_value=True)
    @patch('builtins.open', new_callable=mock_open, read_data=DEFAULT_JSON)
    def test_load(self, mock_open, mock_exists):
        SETTINGS_PATH='/home/user/settings.json'
        s = Settings(SETTINGS_PATH)
        mock_open.assert_called_once_with(SETTINGS_PATH)
        mock_exists.assert_called_once()

        ap = s.aperture
        self.assertEqual(ap.r, 2*u.arcsec)
        self.assertEqual(ap.r_in, 4*u.arcsec)
        self.assertEqual(ap.r_out, 8*u.arcsec)

        self.assertSequenceEqual(s.bands, ['X', 'Y'])

        self.assertEqual(s.get_comp(('X', 'Y')), 'comp-star-1')
        self.assertEqual(s.get_check(('X', 'Y')), 'check-star-2')
        self.assertSetEqual(s.get_rejected(('X', 'Y')), set([1,2,3]))
