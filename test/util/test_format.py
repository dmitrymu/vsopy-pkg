from vsopy.util import default_table_format
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column
import astropy.units as u
import unittest

class DefaultTableFormatterTest(unittest.TestCase):
    def test_simple_formatting(self):
        table = QTable()
        table['flux'] = [1, 2.5, 3.11]
        table['snr'] = [0.1, 1.23, 42.003]
        table['peak'] = [0.01, 0.99999, 0.42]
        table['mag_err'] = Column([(4.2, 0.42), (4.24, 0.4242), (42, 0.4)],
                                  unit=u.mag,
                                  dtype=[('mag', 'f4'), ('err', 'f4')])

        formatted_table = default_table_format(table)

        EXPECTED = '''
flux snr   peak  mag_err [mag, err]
                        mag       \x20
---- ---- ------ ------------------
   1  0.1   1.0%       4.20 ± 0.420
   2  1.2 100.0%       4.24 ± 0.424
   3 42.0  42.0%      42.00 ± 0.400
'''
        self.assertEqual(str(formatted_table), EXPECTED.strip())

    def test_coord_formatting(self):
        table = QTable()
        table['radec2000'] = SkyCoord([[10.25, 20.5] * u.deg, [30.011, 40.0001] * u.deg])

        EXPECTED = '''       radec2000      \x20
        deg,deg       \x20
-----------------------
10°15′00.0″ 20°30′00.0″
30°00′39.6″ 40°00′00.4″'''
        formatted_table = default_table_format(table)
        self.assertEqual(str(formatted_table), EXPECTED)
