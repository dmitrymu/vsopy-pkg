import astropy.units as u
import numpy as np
import unittest

from astropy.table import Column, QTable
from astropy.time import Time
from types import SimpleNamespace
from vsopy.phot import estimate_photometry_quality
from vsopy.util import MagErrDtype


def mock_provider(bands, offsets, residuals):
    star_baselines = np.array([10., 11., 12., 13.])
    auids = [f'star-{i}' for i in range(len(star_baselines))]
    image_ids = np.arange(len(offsets))
    rows = [
        (
            image_index,
            auid,
            baseline + offset + residuals[image_index, star_index],
            .02,
        )
        for image_index, (_, offset) in enumerate(zip(bands, offsets))
        for star_index, (auid, baseline) in enumerate(
            zip(auids, star_baselines)
        )
    ]
    rows.extend([(0, 'target', 9.5, .01), (1, 'target', 9.5, .01)])

    return SimpleNamespace(
        images_=QTable({
            'image_id': image_ids,
            'filter': bands,
            'time': Time('2025-01-01') + image_ids * u.min,
            'airmass': 1 + image_ids * .1,
        }),
        measured_=QTable({
            'image_id': [row[0] for row in rows],
            'auid': [row[1] for row in rows],
            'M': Column(
                [(row[2], row[3]) for row in rows],
                unit=u.mag,
                dtype=MagErrDtype,
            ),
        }),
        target_auid='target',
    )


class PhotometryQualityTest(unittest.TestCase):

    def test_estimates_metrics(self):
        offsets = np.array([-.3, -.1, .1, .3])
        image_pattern = np.array([1., 2., -1., -2.])
        star_pattern = np.array([1., -1., 1., -1.])
        residuals = .01 * np.outer(image_pattern, star_pattern)
        provider = mock_provider(['V'] * 4, offsets, residuals)

        result = estimate_photometry_quality(provider)

        np.testing.assert_allclose(
            result['offset_mag'].value, offsets, atol=3e-7
        )
        np.testing.assert_allclose(
            result['residual_scatter_mag'].value,
            np.sqrt(4 / 3) * .01 * np.abs(image_pattern),
            atol=3e-7,
        )
        np.testing.assert_allclose(
            result['uncertainty_consistency'],
            image_pattern ** 2 / 3,
            atol=2e-5,
        )
        expected_direction = []
        for offset, factor in zip(offsets - offsets[0], image_pattern):
            pre_offset = offset + .01 * factor * star_pattern
            expected_direction.append(
                np.sum(pre_offset)
                / np.sqrt(len(pre_offset) * np.sum(pre_offset ** 2))
            )
        np.testing.assert_allclose(
            result['residual_coherence'], expected_direction, atol=2e-5
        )
        self.assertIn('time', result.colnames)
        self.assertIn('airmass', result.colnames)

    def test_excludes_additional_stars(self):
        provider = mock_provider(
            ['V'] * 4, np.zeros(4), np.zeros((4, 4))
        )

        result = estimate_photometry_quality(
            provider, exclude=['star-3']
        )

        self.assertEqual(len(result), 4)

    def test_single_star_has_no_variance_metrics(self):
        provider = mock_provider(
            ['V'] * 4, np.zeros(4), np.zeros((4, 4))
        )

        result = estimate_photometry_quality(
            provider, exclude=['star-1', 'star-2', 'star-3']
        )

        self.assertTrue(
            np.isnan(result['residual_scatter_mag'].value).all()
        )
        self.assertTrue(np.isnan(result['uncertainty_consistency']).all())


if __name__ == '__main__':
    unittest.main()
