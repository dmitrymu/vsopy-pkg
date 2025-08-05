import astropy.units as u
import numpy as np
from astropy.units import Quantity
from collections.abc import Callable


def ccd_gain(score, full_well, adc_bits):
    scale = (1 << adc_bits) * u.adu
    return full_well / scale / (10 ** (score / (20 * u.db)))


def UNITY_XFM(g, full_well, adc_bits): return 1.0*u.electron/u.adu


def ASI_XFM(g, full_well, adc_bits): return ccd_gain(
    g * u.db / 10, full_well, adc_bits)


class Camera:
    """Camera class representing a CCD camera with specific characteristics.
    """

    def __init__(self, full_well: Quantity, adc_bits: int,
                 xfm: Callable[[float, Quantity, int], Quantity] = UNITY_XFM,
                 read_noise: tuple[list[float], list[Quantity]] = None):
        """Create a Camera instance.

        :param full_well: Full well capacity in electrons
        :type full_well: :py:class:`~astropy.units.Quantity`
        :param adc_bits: Number of bits in the ADC (Analog-to-Digital Converter)
        :type adc_bits: int
        :param xfm: Function to convert gain from camera units to electrons, defaults to unity transformation
        :type xfm: function, optional
        :param read_noise: Read noise characteristics, defaults to None
        :type read_noise: tuple[list[float], list[Quantity]] or None, optional
        """
        self.full_well_ = full_well
        self.adc_bits_ = adc_bits
        self.xfm_ = xfm
        self.read_noise_ = read_noise

    def gain_to_e(self, score: float) -> Quantity:
        """ Convert gain score to electrons.

        :param score: Gain score in camera units
        :type score: float
        :return: Gain in electrons
        :rtype: :py:class:`~astropy.units.Quantity`
        """
        return self.xfm_(score, self.full_well_, self.adc_bits_)

    def read_noise(self, score: float) -> Quantity:
        """Camera read noise in electrons for a given gain score.

        :param score: Gain score in camera units
        :type score: float
        :return: read noise in electrons, or None if read noise curve was not defined
        :rtype: :py:class:`~astropy.units.Quantity`
        """
        return np.interp(score, self.read_noise_[0], self.read_noise_[1]) \
            if self.read_noise_ else None

    @property
    def adu_scale(self) -> int:
        """Scale factor to convert pixel value to ADU count.

        :return: :math:`2^{16 - N}`, where :math:`N` is the number of bits in the ADC.
        :rtype:  int
        """
        return 1 << (16 - self.adc_bits_)

    @property
    def max_adu(self) -> int:
        """Maximum ADU value for the camera.

        :return: :math:`2^{N}`, where :math:`N` is the number of bits in the ADC.
        :rtype:  int
        """
        return (1 << self.adc_bits_) * u.adu


class CameraRegistry:
    """Registry of known cameras

    This class provides access to the registry of known cameras by a camera name.
    The contents of the registry are hardcoded.
    """
    KNOWN_CAMERAS = {
        'ZWO CCD ASI533MM Pro': Camera(50_000 * u.electron,
                                       14,
                                       ASI_XFM,
                                       read_noise=([0, 50, 99, 100, 200, 400],
                                                   [3.8, 3.4, 3.0, 1.46, 1.27, 0.9] * u.electron))
    }

    def __init__(self):
        pass

    @staticmethod
    def get(name) -> Camera:
        """Get a camera by its name from the registry.

        :param name: Name of the camera to retrieve
        :type name: str
        :return: Camera object corresponding to the name, or None if not found
        :rtype: :py:class:`Camera`
        """
        return CameraRegistry.KNOWN_CAMERAS[name] \
            if name in CameraRegistry.KNOWN_CAMERAS \
            else None
