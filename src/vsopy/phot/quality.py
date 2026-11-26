import numpy as np
import statsmodels.api as sm

from astropy import units as u
from astropy.table import Column, QTable, join


def estimate_photometry_quality(provider, exclude=()) -> QTable:
    """Estimate consistency of instrumental photometry within each image.

    The function consumes magnitudes already measured by the photometry
    pipeline. It excludes the session target and fits comparison-star
    measurements separately for every filter using:

    .. math::

        m_{i,s} = \\mu_s + B_i + \\epsilon_{i,s}

    where :math:`\\mu_s` is the instrumental baseline of star :math:`s`,
    :math:`B_i` is the common offset of image :math:`i`, and
    :math:`\\epsilon_{i,s}` is the residual. The fit uses inverse-variance
    weights :math:`1/\\sigma_{i,s}^2`.

    Star baselines and image offsets have no absolute zero point in this
    model. ``offset_mag`` is therefore centered within each filter using a
    measurement-count-weighted session mean.

    :param provider: provider containing ``images_``, ``measured_``, and
                     ``target_auid``.
    :type provider: ~vsopy.phot.BatchDataProvider
    :param exclude: additional comparison-star AUIDs to omit.
    :type exclude: iterable[str]
    :return: one row per image with the original image metadata and:

        * ``offset_mag`` — centered common image offset in magnitudes.
        * ``residual_scatter_mag`` — independent comparison-star scatter:

          .. math::

              \\sqrt{\\frac{\\sum_s \\epsilon_{i,s}^2}{N_i-1}}

        * ``uncertainty_consistency`` — residual scatter normalized by the
          reported uncertainties:

          .. math::

              \\frac{1}{N_i-1}\\sum_s
              \\left(\\frac{\\epsilon_{i,s}}{\\sigma_{i,s}}\\right)^2

          A value near one indicates consistency between scatter and reported
          uncertainty.
        * ``residual_coherence`` — bounded signed indicator from the
          inverse-variance-weighted residuals before removal of the image
          offset. Values near zero indicate no dominant common direction;
          values near either end of ``[-1, 1]`` indicate coherent residuals.

        Scatter and uncertainty-consistency values are ``NaN`` when an image
        has only one valid comparison star because their denominator is zero.
    :rtype: ~astropy.table.QTable
    """
    photometry = join(
        provider.images_['image_id', 'filter'],
        provider.measured_['image_id', 'auid', 'M'],
        keys='image_id',
    )

    excluded = set(exclude)
    if provider.target_auid is not None:
        excluded.add(provider.target_auid)
    if excluded:
        photometry = photometry[
            ~np.isin(photometry['auid'], list(excluded))
        ]

    photometry['M_err'] = photometry['M']['err'].value
    photometry['M'] = photometry['M']['mag'].value

    rows = []
    for band in np.unique(photometry['filter']):
        band_data = photometry[photometry['filter'] == band]
        valid = (
            np.isfinite(band_data['M'])
            & np.isfinite(band_data['M_err'])
            & (band_data['M_err'] > 0)
        )
        band_data = band_data[valid]
        if len(band_data) == 0:
            continue

        star_ids = np.unique(band_data['auid'])
        image_ids = np.unique(band_data['image_id'])
        star_columns = np.column_stack([
            band_data['auid'] == star_id
            for star_id in star_ids[1:]
        ]) if len(star_ids) > 1 else np.empty((len(band_data), 0))
        image_columns = np.column_stack([
            band_data['image_id'] == image_id
            for image_id in image_ids[1:]
        ]) if len(image_ids) > 1 else np.empty((len(band_data), 0))
        design = np.column_stack([
            np.ones(len(band_data)),
            star_columns,
            image_columns,
        ]).astype(float)

        result = sm.WLS(
            band_data['M'],
            design,
            weights=1 / band_data['M_err'] ** 2,
        ).fit()
        band_data['residual'] = result.resid

        offset_count = len(image_ids) - 1
        offset_coefficients = (
            result.params[-offset_count:]
            if offset_count else []
        )
        raw_offsets = np.concatenate(([0.0], offset_coefficients))
        image_offsets = dict(zip(image_ids, raw_offsets))

        image_counts = np.array([
            np.count_nonzero(band_data['image_id'] == image_id)
            for image_id in image_ids
        ])
        offset_center = np.average(
            np.fromiter(image_offsets.values(), dtype=float),
            weights=image_counts,
        )

        for image_id in image_ids:
            image_data = band_data[band_data['image_id'] == image_id]
            residual = image_data['residual'].data
            errors = image_data['M_err'].data
            weights = 1 / errors ** 2
            pre_offset_residual = residual + image_offsets[image_id]
            denominator = np.sqrt(
                np.sum(weights)
                * np.sum(weights * pre_offset_residual ** 2)
            )
            count = len(image_data)

            rows.append({
                'image_id': image_id,
                'filter': band,
                'offset_mag': image_offsets[image_id] - offset_center,
                'residual_scatter_mag': (
                    np.sqrt(np.sum(residual ** 2) / (count - 1))
                    if count > 1 else np.nan
                ),
                'uncertainty_consistency': (
                    np.sum((residual / errors) ** 2) / (count - 1)
                    if count > 1 else np.nan
                ),
                'residual_coherence': (
                    np.sum(weights * pre_offset_residual) / denominator
                    if denominator > 0 else np.nan
                ),
            })

    if not rows:
        return QTable()

    quality = QTable(rows)
    for name in ('offset_mag', 'residual_scatter_mag'):
        quality[name] = Column(quality[name], unit=u.mag)

    return join(
        quality,
        provider.images_,
        keys=['image_id', 'filter'],
    )
