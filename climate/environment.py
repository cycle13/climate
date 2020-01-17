import numpy as np

from .constants import KELVIN

def sky_emissivity(dew_point_temperature, total_sky_cover):
    """
    Calculate the sky emissivity using the Clark-Allen method used by EnergyPlus.

    Clark and Allen report that the standard estimation error for the atmospheric long-wave radiation was 10 𝑊𝑊/𝑚𝑚2.

    Parameters
    ----------
    dew_point_temperature : float
        Dew-point temperature of air in degrees C
    total_sky_cover : int
        Total sky cover in tenths. Conversion to decimal 0-1 happens within the method.

    Returns
    -------
    sky_emissivity : float

    """

    sky_emissivity = (0.787 + 0.764 * np.log((dew_point_temperature + KELVIN) / KELVIN)) * (
            1 + 0.0224 * total_sky_cover - 0.0035 * np.power(total_sky_cover, 2) + 0.0028 * np.power(
        total_sky_cover, 3))
    return sky_emissivity