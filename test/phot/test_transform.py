import astropy.units as u
import numpy as np
import unittest

from astropy.table import QTable, Column
from vsopy.mock import mock_measure
from vsopy.phot import apply_simple_transform, create_simple_transform, SimpleTransform
from vsopy.util import  MagErrDtype, ValErr

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

    def test_simple_transform(self):
        np.random.seed(42)  # for reproducibility
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

        A_c, B_c = (ValErr(10.1, 0.03), ValErr(9.7, 0.03))
        a_c, b_c = (ValErr(7.40, 0.061), ValErr(7.05, 0.055))
        a_t, b_t = (ValErr(8.50, 0.47), ValErr(8.77, 0.052))

        A_t1, B_t2, B_t1, A_t2 = apply_simple_transform(
            xfm, A_c, B_c, a_c, b_c, a_t, b_t
        )

        self.assertAlmostEqual(A_t1.mag, 10.6, places=1)
        self.assertAlmostEqual(B_t2.mag, 12.1, places=1)
        self.assertAlmostEqual(A_t2.mag, 10.75, places=1)
        self.assertAlmostEqual(B_t1.mag, 11.9, places=1)