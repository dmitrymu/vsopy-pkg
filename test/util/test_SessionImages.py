import unittest
import numpy as np
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u
from numpy.testing import assert_array_equal, assert_array_almost_equal
from pathlib import Path
from vsopy.util.SessionImages import session_image_list, batch_session_images
from unittest.mock import patch, MagicMock, PropertyMock

def mock_dir(name):
    """Helper function to create a mock directory."""
    result = MagicMock(spec=Path)
    result.is_dir.return_value = True
    result.name = name
    result.__str__.return_value = name
    result.__truediv__.side_effect = lambda x: mock_dir(str(Path(name) / x))
    return result

class SessionImageListTest(unittest.TestCase):

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

class BatchSessionImagesTest(unittest.TestCase):

    @patch(f"vsopy.util.SessionImages.QTable.read")
    def test_batch_session_images(self, mock_read):
        mock_read.return_value = QTable({
            'image_id': [1, 2, 3, 4, 5, 6],
            'filter': ['V', 'R', 'Ic', 'V', 'R', 'Ic'],
            'time': Time(['2023-10-01T00:00:00', '2023-10-01T00:00:11',
                          '2023-10-01T00:00:34', '2023-10-01T00:01:10',
                          '2023-10-01T00:01:21', '2023-10-01T00:01:44']),
            'exposure': [10, 20, 30, 10, 20, 30] * u.second,
            'airmass': [1.2, 1.201, 1.202, 1.21, 1.211, 1.212],
            'temperature': [-10, -10.1, -9.9, -10.1, -10, -9.9] * u.deg_C,
            'file': ['image1.fits', 'image2.fits', 'image3.fits', 'image4.fits',
                     'image5.fits', 'image6.fits'],
        })


        batches, batch_images = batch_session_images('dummy_path')
        self.assertEqual(len(batches), 2)
        self.assertEqual(len(batch_images), 6)

        expected_batches = QTable({
            'batch_id': [1, 2],
            'time': Time(['2023-10-01T00:00:15', '2023-10-01T00:01:25'],
                          scale='utc', format='isot'),
            'time_range': [34.0, 34.0] * u.second,
            'temperature': [-10.0, -10.0] * u.deg_C,
            'temperature_range': [.2, .2] * u.deg_C,
            'airmass': [1.201, 1.211],
            'airmass_range': [0.002, 0.002]
        })

        assert_array_equal(batches['batch_id'], expected_batches['batch_id'])
        assert_array_almost_equal(batches['time'].jd, expected_batches['time'].jd)
        assert_array_almost_equal(batches['time_range'].value, expected_batches['time_range'].value)
        assert_array_almost_equal(batches['temperature'].value, expected_batches['temperature'].value)
        assert_array_almost_equal(batches['temperature_range'].value, expected_batches['temperature_range'].value)
        assert_array_almost_equal(batches['airmass'].value, expected_batches['airmass'].value)
        assert_array_almost_equal(batches['airmass_range'].value, expected_batches['airmass_range'].value)

        expected_batch_images = QTable({
            'batch_id': [1, 1, 1, 2, 2, 2],
            'image_id': [1, 2, 3, 4, 5, 6]
        })

        assert_array_equal(batch_images['batch_id'], expected_batch_images['batch_id'])
        assert_array_equal(batch_images['image_id'], expected_batch_images['image_id'])