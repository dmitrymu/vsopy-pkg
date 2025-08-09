import astropy.units as u
import numpy as np
import unittest

from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column
from astropy.time import Time
from numpy.testing import assert_array_equal
from pathlib import Path
from vsopy.phot import BatchDataProvider
from vsopy.util import MagErrDtype
from unittest.mock import patch, MagicMock

BATCHES = QTable(dict(
    batch_id=[1, 2],
    time=Time(['2023-01-01T01:00:00', '2023-01-02T01:05:00']),
    temperature=[20.0, 21.0] * u.deg_C,
    airmass=[1.05, 1.25],
    time_range=[120, 115] * u.second,
    temperature_range=[.1, .2] * u.deg_C,
    airmass_range=[.01, .02],
))

BATCHES_IMAGES = QTable(dict(
    batch_id=[1, 1, 2, 2],
    image_id=[1, 2, 3, 4],
))

IMAGES = QTable(dict(
    image_id=[1, 2, 3, 4],
    filter=['A', 'A', 'B', 'B'],
    time=Time(['2023-01-01T01:00:10', '2023-01-01T01:02:20',
               '2023-01-02T01:05:10', '2023-01-02T01:07:20']),
    exposure=[60, 60, 60, 60] * u.second,
    airmass=[1.0, 1.1, 1.2, 1.3],
    temperature=[19.05, 20.05, 20.9, 21.1] * u.deg_C,
    gain=[1.0, 1.0, 1.0, 1.0],
))

MEASURED = QTable(dict(
    image_id=[1, 1, 2, 2, 3, 3, 4, 4],
    auid=['star1', 'star2', 'star1', 'star2', 'star1', 'star2', 'star1', 'star2'],
    M=Column([(10.0, .05), (11.0, .02), (10.5, .03), (11.5, .01),
              (12.0, .04), (13.0, .02), (12.5, .03), (13.5, .01)] * u.mag,
             dtype=MagErrDtype),
    flux=[1000, 800, 900, 700, 600, 500, 400, 300] * (u.electron/u.second),
    snr=[50, 40, 45, 35, 30, 25, 20, 15] * u.db,
    peak=[1200, 1000, 1100, 900, 800, 700, 600, 500],
))

SEQUENCE = QTable(dict(
    auid=['star1', 'star2', 'star3', 'star1', 'star2', 'star3'],
    band=['A', 'A', 'A', 'B', 'B', 'B'],
    M=Column([(10.0, .05), (11.0, .02), (3.0, .02),
              (12.0, .04), (13.0, .02), (4.0, .02)] * u.mag,
             dtype=MagErrDtype),
    ),
    meta={'auid': 'star1'}
)

def mock_session(root):
    session = MagicMock()
    session.batches_file_path.return_value = Path(root) / 'batches.ecsv'
    session.batch_images_file_path.return_value = Path(root) / 'batch_images.ecsv'
    session.images_file_path.return_value = Path(root) / 'images.ecsv'
    session.measured_file_path.return_value = Path(root) / 'measured.ecsv'
    session.sequence_file_path.return_value = Path(root) / 'sequence.ecsv'
    return session

class MeasurePhotometryTest(unittest.TestCase):

    @patch(f"vsopy.phot.batch_data_provider.QTable.read")
    def test_provider(self, mock_read):
        mock_read.side_effect = [BATCHES, BATCHES_IMAGES, IMAGES, MEASURED, SEQUENCE]
        session = mock_session('/home/test')
        provider = BatchDataProvider(session)

        self.assertEqual(provider.target_auid, 'star1')
        sbp = provider.sequence_band_pair(('A', 'B'))
        self.assertSequenceEqual(sbp.colnames, ['auid', 'A', 'B'])
        # star3 should not be included in the sequence
        self.assertEqual(len(sbp), 2)
        assert_array_equal(sbp['auid'], ['star1', 'star2'])