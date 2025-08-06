import astropy.units as u
import numpy as np
import unittest

from astropy.table import QTable, Column
from vsopy.mock import mock_measure
from vsopy.phot import create_simple_transform, SimpleTransform
from vsopy.util import  MagErrDtype

SEQUENCE = QTable(dict(
    auid=[f'star-{i}' for i in range(10)],
    A=Column([(10.0, .04), (10.2, .03), (10.4, .06), (10.6, .06), (10.8, .06),
              (11.0, .04), (11.2, .03), (11.4, .06), (11.6, .06), (11.8, .06)],
             unit=u.mag, dtype=MagErrDtype),
    B=Column([(9.0, .04), (9.4, .03), (9.8, .06), (10.2, .06), (10.6, .06),
              (11.0, .04), (11.4, .03), (11.8, .06), (12.2, .06), (12.6, .06)],
             unit=u.mag, dtype=MagErrDtype)
))

class TransformTest(unittest.TestCase):

    def test_create_simple_transform(self):
        np.random.seed(42)
        measured = mock_measure(SEQUENCE, ('A', 'B'),
                                (0.9, 0.02, 0.05), (1.1, 5.0, 0.05))
        xfm = create_simple_transform(
            SEQUENCE['A']['mag'], SEQUENCE['B']['mag'],
            measured['A'], measured['B']
        )
        self.assertIsInstance(xfm, SimpleTransform)
        self.assertAlmostEqual(xfm.Ta.val, 1.1, places=1)
        #self.assertAlmostEqual(xfm.Ta.err, 0.05, places=2)
        self.assertAlmostEqual(xfm.Tb.val, 1.2, places=1)
        #@self.assertAlmostEqual(xfm.Tb.err, 0.02, places=2)
        self.assertAlmostEqual(xfm.Tab.val, 0.9, places=1)
        #self.assertAlmostEqual(xfm.Tab.err, 0.05/0.9**2, places=2)