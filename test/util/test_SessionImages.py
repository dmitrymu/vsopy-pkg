import unittest
import numpy as np
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u
from numpy.testing import assert_array_equal
from pathlib import Path
from vsopy.util.SessionImages import session_image_list
from unittest.mock import patch, MagicMock, PropertyMock

def mock_dir(name):
    """Helper function to create a mock directory."""
    result = MagicMock(spec=Path)
    result.is_dir.return_value = True
    result.name = name
    result.__str__.return_value = name
    result.__truediv__.side_effect = lambda x: mock_dir(str(Path(name) / x))
    return result

class DefaultImageListTest(unittest.TestCase):

    @patch(f"vsopy.util.SessionImages.ccdp.ImageFileCollection.__init__", return_value=None)
    @patch(f"vsopy.util.SessionImages.ccdp.ImageFileCollection.summary", new_callable=PropertyMock)
    def test_session_image_list(self, mock_summary, mock_image_collection):
        mock_summary.side_effect = [
            QTable({
                'file': ['image_V1.fits', 'image_V2.fits'],
                'date-obs': ['2023-10-01T00:00:00', '2023-10-01T01:00:00'],
                'exptime': [10.0, 20.0],
                'ccd-temp': [-10.0, -12.0],
                'filter': ['V', 'V'],
                'airmass': [1.2, 1.5]
            }),
            QTable({
                'file': ['image_R1.fits', 'image_R2.fits'],
                'date-obs': ['2023-10-01T00:00:10', '2023-10-01T01:00:10'],
                'exptime': [10.0, 20.0],
                'ccd-temp': [-10.0, -12.0],
                'filter': ['R', 'R'],
                'airmass': [1.2, 1.5]
            }),
        ]

        session = MagicMock()
        session.tag = "20231001"
        session.name = "TestSession"
        session.lights_dir = MagicMock()
        session.lights_dir.iterdir.return_value = [
            mock_dir("lights/V"),
            mock_dir("lights/R"),
        ]

        images = session_image_list(session)

        self.assertEqual(len(images), 4)
        self.assertIn('image_id', images.colnames)
        self.assertIn('filter', images.colnames)
        self.assertIn('time', images.colnames)
        self.assertIn('exposure', images.colnames)
        self.assertIn('airmass', images.colnames)
        self.assertIn('temperature', images.colnames)
        self.assertIn('path', images.colnames)
        assert_array_equal(images['filter'], np.array(['V', 'V', 'R', 'R']))
        assert_array_equal(images['image_id'], np.array([1, 2, 3, 4]))
        assert_array_equal(images['time'], Time(['2023-10-01T00:00:00', '2023-10-01T01:00:00',
                                               '2023-10-01T00:00:10', '2023-10-01T01:00:10']))
        assert_array_equal(images['exposure'], np.array([10.0, 20.0, 10.0, 20.0]) * u.second)
        assert_array_equal(images['airmass'], np.array  ([1.2, 1.5, 1.2, 1.5]))
        assert_array_equal(images['temperature'], np.array([-10.0, -12.0, -10.0, -12.0]) * u.deg_C)
        assert_array_equal(images['path'], np.array([
            'lights/V/image_V1.fits',
            'lights/V/image_V2.fits',
            'lights/R/image_R1.fits',
            'lights/R/image_R2.fits'
            ]))
