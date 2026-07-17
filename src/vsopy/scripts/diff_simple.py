import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..')))

import numpy as np
import argparse
import json
from vsopy import phot, data, util
from pathlib import Path
from astropy.table import QTable, join
from collections import namedtuple
from operator import attrgetter
from dataclasses import dataclass

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate differential photometry report for the session')
    parser.add_argument('-O', '--object', type=str,
                        required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str,
                        required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str,
                        required=True, help='Work directory')
    parser.add_argument('--observer', type=str,
                        required=True, help='AAVSO observer code')
    parser.add_argument('--sigma', type=float, default=3.0,
                        help='Number of stdev to filter outliers by check star')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')
    parser.add_argument('--no-lsq', action='store_true',
                        default=False, help='Don''t use least square solver')

    return parser.parse_args()

@dataclass
class BandResult:
    data: QTable
    err: float
    comp: str
    check: str

@dataclass
class BandPairResult:
    a: BandResult
    b: BandResult

def columns_not_m(table, m):
    cnames = table.column_names
    exclude = [m, f'check {m}']
    return table[[name for name in cnames if name not in exclude]]

def measure_band_pair(provider, band, settings):
        comp_auid = settings.get_comp(band)
        check_auid = settings.get_check(band)
        targ = phot.batch_diff_photometry(
            provider,
            band,
            comp_auid,
            provider.target_auid
        )
        check = phot.batch_diff_photometry(
            provider,
            band,
            comp_auid,
            check_auid
        )
        check.rename_columns([band[0], band[1]], [f'check {band[0]}', f'check {band[1]}'])
        ab = join(join(targ, provider.batches_, keys='batch_id'), check, keys='batch_id')
        a_err = np.linalg.norm(ab[band[0]]['err'])
        b_err = np.linalg.norm(ab[band[1]]['err'])
        return BandPairResult(
            a = BandResult(columns_not_m(band[1]), a_err, comp_auid, check_auid),
            b = BandResult(columns_not_m(band[0]), b_err, comp_auid, check_auid)
        )


def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    layout = util.WorkLayout(args.work_dir)
    session_layout = layout.get_session(session)
    settings = util.Settings(session_layout.settings_file_path)

    provider = phot.BatchDataProvider(session_layout)
    comp_star = provider.get_star(settings.get_comp(band))
    check_star = provider.get_star(settings.get_check(band))

    def trim_time(data, begin, end):
        flt = data['']

    def reject_outliers(data, band, check_mag, sigma=3.0):
        dc = data[f'check {band}']['mag'] - check_mag
        dc_mean = np.mean(dc)
        dc_std = np.std(dc)
        flt = np.abs(dc - dc_mean) < sigma*dc_std
        return data[flt], np.count_nonzero(~flt)

    bands = util.band_pairs(util.ordered_bands(settings.bands))
    result = {}
    for band in bands:

        ab = measure_band_pair(provider, band, settings)

        result.setdefault(band[0], []).append(ab.a)
        result.setdefault(band[1], []).append(ab.b)


    with open(session_layout.root_dir / 'report-simple.txt', mode='w') as f:
      report = data.AavsoReport(f,
                                provider.target_auid,
                                provider.chart_id,
                                args.observer)
      report.header()
      for band in settings.bands:
        comp = comp_star[comp_star['band'] == band]
        check = check_star[comp_star['band'] == band]
        item = min(result[band], key=attrgetter('err'))
        dt = item.data
        print(f"Selected {band} data with total error {item.err}")
        passed, rejected_count = reject_outliers(dt, band, comp['M']['mag'][0], args.sigma)
        print(f"Rejected {rejected_count} data points outside {args.sigma} sigma")
        report.body(passed, band, comp, check)

    return 0

# Example: python3 diff_simple.py -O RR_Lyr -t 20230704 -w /home/user/work --observer XYZ


if __name__ == '__main__':
    sys.exit(main())
