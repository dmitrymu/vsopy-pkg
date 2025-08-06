import astropy.units as u # type: ignore
import numpy as np

from ..util import MagErr, ValErr, MagErrDtype, ValErrDtype
from astropy.table import QTable, Column, join # type: ignore
from scipy import stats as sst
from typing import NamedTuple

class SimpleTransform(NamedTuple):
    """Simple linear transformation coefficients for two bands.
    """
    Ta: ValErr
    """Transformation coefficient for band A"""
    Tb: ValErr
    """Transformation coefficient for band B"""
    Tab: ValErr
    """Transformation coefficient for color index (A-B)"""

def create_simple_transform(A:np.ndarray, B:np.ndarray,
                            a:np.ndarray, b:np.ndarray) -> SimpleTransform:
    """Create transform from star photometry for two images

    Given an ensemble of stars with both standard (:math:`A`, :math:`B`) and instrumental
    (:math:`a`, :math:`b`) magnitudes known for two bands, we fit two linear regressions:

    .. math::
        :label: create.1

        a - b = T_{ab} (A - B) + Z_{ab}

    .. math::
        :label: create.2

        A - a = T_a (A - B) + Z_a

    Instead of :eq:`create.2`, we can also use the following equation:

    .. math::
        :label: create.3

        B - b = T_b (A - B) + Z_b

    :math:`T_{ab}` determines transformation of instrumental color index to standard
    color index.  :math:`T_a` corrects transformation from instrumental to standard
    magnitude with respect to color index. See :py:func:`apply_simple_transform`
    for details.
    :math:`T_b` is similar to :math:`T_a` but for band B.

    The regression coefficients are calculated using :py:func:`scipy.stats.linregress`.

    :param A: standard magnitudes in band A
    :param B: standard magnitudes in band B
    :param a: instrumental magnitudes in band A
    :param b: instrumental magnitudes in band B
    :type A: :py:class:`~numpy.ndarray`
    :type B: :py:class:`~numpy.ndarray`
    :type a: :py:class:`~numpy.ndarray`
    :type b: :py:class:`~numpy.ndarray`

    :return: transform coefficients :py:class:`SimpleTransform`
    :rtype: SimpleTransform
    :raises Exception: if any of the regression slopes has zero standard error

    See `Gary B.L., CCD Transformation Equations for Use with Single-Image Photometry
    <https://researchers.usask.ca/gordon-sarty/documents/astronomy/transf.pdf>`__
    for details.
    """
    AB = A - B
    ab = a - b
    Aa =A - a
    Bb = B - b
    reg_ab = sst.linregress(AB, ab)
    reg_Aa = sst.linregress(AB, Aa)
    reg_Bb = sst.linregress(AB, Bb)

    if reg_Aa.stderr == 0 or reg_Bb.stderr == 0 or reg_ab.stderr == 0:
        raise Exception("Zero stderr")

    result = SimpleTransform(
        ValErr(reg_Aa.slope, reg_Aa.stderr),
        ValErr(reg_Bb.slope, reg_Bb.stderr),
        ValErr(1/reg_ab.slope, reg_ab.stderr))
    return result

def apply_simple_transform(xfm:SimpleTransform,
                           A_c:MagErr, B_c:MagErr,
                           a_c:MagErr, b_c:MagErr,
                           a_t:MagErr, b_t:MagErr) -> tuple[MagErr, MagErr, MagErr, MagErr]:
    """Simple differential photometry

    Differential photometry is a method of measuring the magnitude of a
    target star with respect to a comparison star, for which instrumental
    magnitude is known. The simple method assumes that the we have a couple of
    images in different spectral bands, and each image contains several stars
    with known standard magnitudes. We select one star as a comparison star,
    another as a check star, and the rest for calculating the transform
    using the method described in :py:func:`create_simple_transform`.

    Let we have images in bands A and B, and the following magnitudes are known:

    * :math:`A_c` and :math:`B_c` - standard magnitudes of the comparison star;
    * :math:`a_c` and :math:`b_c` - instrumental magnitudes of the comparison star;
    * :math:`a_t` and :math:`b_t` - instrumental magnitudes of the target star.

    Let the transform coefficients are :math:`T_a`, :math:`T_b`, and :math:`T_{ab}`.

    The standard magnitudes of the target star are calculated as follows:

    .. math::
        :label: apply.1

        C_t = A_t - B_t = (A_c - B_c) + T_{ab} ((a_t-b_t) - (a_c-b_c))

    .. math::
        :label: apply.2

        A_{t1} = a_t + (A_c-a_c) + T_a C_t

    .. math::
        :label: apply.3

        B_{t1} = A_t - C_t

    We can also calculate transformed magnitudes in the opposite order:

    .. math::
        :label: apply.4

        B_{t2} = b_t + (B_c-b_c) + T_b C_t`


    .. math::
        :label: apply.5

        A_{t2} = B_{t2} + C_t

    Because of chain nature of the equations, the magnitude which is calculated
    first has less uncertainty than the second one. This method returns
    all four target magnitudes, in the order
    :math:`A_{t1}, B_{t2}, A_{t2}, B_{t1}`.

    :param xfm: Transform coefficients
    :type xfm: SimpleTransform
    :param A_c: Standard magnitude of the comparison star in band A
    :type A_c: MagErr
    :param B_c: Standard magnitude of the comparison star in band B
    :type B_c: MagErr
    :param a_c: Instrumental magnitude of the comparison star in band A
    :type a_c: MagErr
    :param b_c: Instrumental magnitude of the comparison star in band B
    :type b_c: MagErr
    :param a_t: Instrumental magnitude of the target star in band A
    :type a_t: MagErr
    :param b_t: Instrumental magnitude of the target star in band B
    :type b_t: MagErr
    :return: Transformed standard magnitudes of the target star in bands A and B,
             in the order :math:`A_{t1}, B_{t2}, A_{t2}, B_{t1}`.
    :rtype: tuple[MagErr, MagErr, MagErr, MagErr]

    See `Gary B.L., CCD Transformation Equations for Use with Single-Image Photometry
    <https://researchers.usask.ca/gordon-sarty/documents/astronomy/transf.pdf>`__
    for details.
    """
    def transform(T_a, T_ab, A_c, B_c, a_c, b_c, a_t, b_t):
        Ta, Ta_err = T_a
        Tab, Tab_err = T_ab
        Ac, Ac_err = A_c
        Bc, Bc_err = B_c
        ac, ac_err = a_c
        bc, bc_err = b_c
        at, at_err = a_t
        bt, bt_err = b_t

        # c_t = a_t - b_t
        atbt = at-bt
        # c_C = a_c - b_c
        acbc = ac-bc

        atbtacbc = atbt - acbc
        atbtacbc_err = np.sqrt(at_err**2 + bt_err**2 + ac_err**2 + bc_err**2)
        Tab_abab = Tab * (atbt - acbc)
        Tab_abab_err = Tab * atbtacbc * \
            np.sqrt((Tab_err/Tab)**2 + (atbtacbc_err/atbtacbc)**2)
        # (1)
        AtBt = (Ac-Bc) + Tab * (atbt - acbc)
        AtBt_err = np.sqrt(Ac_err**2 + Bc_err**2 + Tab_abab_err**2)
        # (2)
        Acac = Ac - ac
        Ta_Tab_err = Ta * Tab_abab * \
            np.sqrt((Ta_err/Ta)**2 + (Tab_abab_err/Tab_abab)**2)
        At = at + Acac + Ta * Tab * (atbt - acbc)
        At_err = np.sqrt(at_err**2 + Ac_err**2 + ac_err**2 + Ta_Tab_err**2)

        # (3)
        Bt = At - AtBt
        Bt_err = np.sqrt(At_err**2 + AtBt_err**2)

        return MagErr(At, At_err), MagErr(Bt, Bt_err)

    Ta, Tb, Tab = xfm
    r1 =  transform(Ta, Tab, A_c, B_c, a_c, b_c, a_t, b_t)
    r2 =  transform(Tb, Tab, B_c, A_c, b_c, a_c, b_t, a_t)
    return r1[0], r2[0], r2[1], r1[1]

def batch_create_simple_transform(provider, bands):

    data = provider.batch_and_sequence_band_pair(bands)
    grouped = data.group_by(['batch_id'])

    xfm = [create_simple_transform(batch[bands[0]]['mag'], batch[bands[1]]['mag'],
                                   batch[provider.instr(bands[0])]['mag'], batch[provider.instr(bands[1])]['mag'])
            for batch in grouped.groups]
    xfm_dtype = [('Ta', ValErrDtype), ('Tb', ValErrDtype), ('Tab', ValErrDtype)]
    return QTable({
        'batch_id': grouped.groups.keys['batch_id'],
        'xfm': Column(xfm, dtype = xfm_dtype)
    })

def batch_apply_simple_transform(provider, xfm, bands, comparison_auid, target_auid=None):

    comp = provider.batch_comp_star(bands, comparison_auid)
    targ = provider.batch_target_star(bands, target_auid if target_auid is not None else provider.target_auid)
    xfm_input = join(join(comp, targ, 'batch_id'), xfm, 'batch_id')
    transformed = [apply_simple_transform(row['xfm'], row[bands[0]].value, row[bands[1]].value,
                            row[provider.instr(bands[0])].value, row[provider.instr(bands[1])].value,
                            row[provider.targ(bands[0])].value, row[provider.targ(bands[1])].value
                            ) for row in xfm_input]
    a, b, _, _ =zip(*transformed)
    return QTable({
        'batch_id': xfm_input['batch_id'],
        bands[0]: Column(list(a), dtype = MagErrDtype, unit=u.mag),
        bands[1]: Column(list(b), dtype = MagErrDtype, unit=u.mag)
        })

def batch_diff_photometry(provider, bands, comparison_auid, target_auid=None):
    """Magnitude transformation for differential photometry.

    It uses simplified approach which ignores second order extinction.
    See Gary B.L., CCD Transformation Equations for Use with Single-Image Photometry.
    This transform must be calculated for each pair of images separately and applied
    to the instrumental magnitudes extracted from those images.
    """

    xfm = batch_create_simple_transform(provider, bands)
    return batch_apply_simple_transform(provider, xfm, bands, comparison_auid, target_auid)

