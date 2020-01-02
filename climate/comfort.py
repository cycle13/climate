# TODO: Add comfort metrics into the mix here!
from climate.constants import KELVIN

from psychrolib import GetSatVapPres

import warnings
import numpy as np
import pandas as pd

def universal_thermal_climate_index(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity):
    """
    Return the approximate Universal Thermal Climate Index for the input criteria.
    :param air_temperature: Air temperature (C)
    :param mean_radiant_temperature: Mean radiant temperature (C)
    :param air_velocity: Air velocity (m/s)
    :param relative_humidity: Relative humidity
    :return utci_approx: Approximate UTCI temperature (C)
    """

    # Convert inputs to numpy arrays for matricizations!
    air_temperature = np.array(air_temperature) if isinstance(air_temperature, list) else np.array([air_temperature])
    mean_radiant_temperature = np.array(mean_radiant_temperature) if isinstance(mean_radiant_temperature, list) else np.array([mean_radiant_temperature])
    air_velocity = np.array(air_velocity) if isinstance(air_velocity, list) else np.array([air_velocity])
    relative_humidity = np.array(relative_humidity) if isinstance(relative_humidity, list) else np.array([relative_humidity])

    # Check for valid application limits and raise warning if these are exceeded
    if np.any(air_temperature < -50) | np.any(air_temperature > 50):
        warnings.warn('Air temperature out of range. -50C and 50C in order for this method to be accurate.', Warning)
    if np.any(mean_radiant_temperature < 30) | np.any(mean_radiant_temperature > 70):
        warnings.warn('Mean radiant temperature should be between 30C and 70C in order for this method to be accurate.', Warning)
    if np.any(air_velocity < 0.5) | np.any(air_velocity > 17):
        warnings.warn('Air velocity should be between 0.5m/s and 17m/s in order for this method to be accurate. Values outside this range will be replaced by the nearest valid value', Warning)
    air_velocity = np.where(air_velocity < 0.5, 0.5, np.where(air_velocity > 17, 17, air_velocity))

    saturated_vapor_pressure = np.apply_along_axis(GetSatVapPres, 0, air_temperature) * 0.01 * (relative_humidity / 100.0) / 10.0
    mean_radiant_temperature_delta = mean_radiant_temperature - air_temperature

    utci_approx = air_temperature + \
                  0.607562052 + \
                  -0.0227712343 * air_temperature + \
                  8.06470249E-4 * np.power(air_temperature, 2) + \
                  -1.54271372E-4 * np.power(air_temperature, 3) + \
                  -3.24651735E-6 * np.power(air_temperature, 4) + \
                  7.32602852E-8 * np.power(air_temperature, 5) + \
                  1.35959073E-9 * np.power(air_temperature, 6) + \
                  -2.25836520 * air_velocity + \
                  0.0880326035 * air_temperature * air_velocity + \
                  0.00216844454 * np.power(air_temperature, 2) * air_velocity + \
                  -1.53347087E-5 * np.power(air_temperature, 3) * air_velocity + \
                  -5.72983704E-7 * np.power(air_temperature, 4) * air_velocity + \
                  -2.55090145E-9 * np.power(air_temperature, 5) * air_velocity + \
                  -0.751269505 * np.power(air_velocity, 2) + \
                  -0.00408350271 * air_temperature * np.power(air_velocity, 2) + \
                  -5.21670675E-5 * np.power(air_temperature, 2) * np.power(air_velocity, 2) + \
                  1.94544667E-6 * np.power(air_temperature, 3) * np.power(air_velocity, 2) + \
                  1.14099531E-8 * np.power(air_temperature, 4) * np.power(air_velocity, 2) + \
                  0.158137256 * np.power(air_velocity, 3) + \
                  -6.57263143E-5 * air_temperature * np.power(air_velocity, 3) + \
                  2.22697524E-7 * np.power(air_temperature, 2) * np.power(air_velocity, 3) + \
                  -4.16117031E-8 * np.power(air_temperature, 3) * np.power(air_velocity, 3) + \
                  -0.0127762753 * np.power(air_velocity, 4) + \
                  9.66891875E-6 * air_temperature * np.power(air_velocity, 4) + \
                  2.52785852E-9 * np.power(air_temperature, 2) * np.power(air_velocity, 4) + \
                  4.56306672E-4 * np.power(air_velocity, 5) + \
                  -1.74202546E-7 * air_temperature * np.power(air_velocity, 5) + \
                  -5.91491269E-6 * np.power(air_velocity, 6) + \
                  0.398374029 * mean_radiant_temperature_delta + \
                  1.83945314E-4 * air_temperature * mean_radiant_temperature_delta + \
                  -1.73754510E-4 * np.power(air_temperature, 2) * mean_radiant_temperature_delta + \
                  -7.60781159E-7 * np.power(air_temperature, 3) * mean_radiant_temperature_delta + \
                  3.77830287E-8 * np.power(air_temperature, 4) * mean_radiant_temperature_delta + \
                  5.43079673E-10 * np.power(air_temperature, 5) * mean_radiant_temperature_delta + \
                  -0.0200518269 * air_velocity * mean_radiant_temperature_delta + \
                  8.92859837E-4 * air_temperature * air_velocity * mean_radiant_temperature_delta + \
                  3.45433048E-6 * np.power(air_temperature, 2) * air_velocity * mean_radiant_temperature_delta + \
                  -3.77925774E-7 * np.power(air_temperature, 3) * air_velocity * mean_radiant_temperature_delta + \
                  -1.69699377E-9 * np.power(air_temperature, 4) * air_velocity * mean_radiant_temperature_delta + \
                  1.69992415E-4 * np.power(air_velocity, 2) * mean_radiant_temperature_delta + \
                  -4.99204314E-5 * air_temperature * np.power(air_velocity, 2) * mean_radiant_temperature_delta + \
                  2.47417178E-7 * np.power(air_temperature, 2) * np.power(air_velocity, 2) * mean_radiant_temperature_delta + \
                  1.07596466E-8 * np.power(air_temperature, 3) * np.power(air_velocity, 2) * mean_radiant_temperature_delta + \
                  8.49242932E-5 * np.power(air_velocity, 3) * mean_radiant_temperature_delta + \
                  1.35191328E-6 * air_temperature * np.power(air_velocity, 3) * mean_radiant_temperature_delta + \
                  -6.21531254E-9 * np.power(air_temperature, 2) * np.power(air_velocity, 3) * mean_radiant_temperature_delta + \
                  -4.99410301E-6 * np.power(air_velocity, 4) * mean_radiant_temperature_delta + \
                  -1.89489258E-8 * air_temperature * np.power(air_velocity, 4) * mean_radiant_temperature_delta + \
                  8.15300114E-8 * np.power(air_velocity, 5) * mean_radiant_temperature_delta + \
                  7.55043090E-4 * np.power(mean_radiant_temperature_delta, 2) + \
                  -5.65095215E-5 * air_temperature * np.power(mean_radiant_temperature_delta, 2) + \
                  -4.52166564E-7 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 2) + \
                  2.46688878E-8 * np.power(air_temperature, 3) * np.power(mean_radiant_temperature_delta, 2) + \
                  2.42674348E-10 * np.power(air_temperature, 4) * np.power(mean_radiant_temperature_delta, 2) + \
                  1.54547250E-4 * air_velocity * np.power(mean_radiant_temperature_delta, 2) + \
                  5.24110970E-6 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 2) + \
                  -8.75874982E-8 * np.power(air_temperature, 2) * air_velocity * np.power(mean_radiant_temperature_delta, 2) + \
                  -1.50743064E-9 * np.power(air_temperature, 3) * air_velocity * np.power(mean_radiant_temperature_delta, 2) + \
                  -1.56236307E-5 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) + \
                  -1.33895614E-7 * air_temperature * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) + \
                  2.49709824E-9 * np.power(air_temperature, 2) * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) + \
                  6.51711721E-7 * np.power(air_velocity, 3) * np.power(mean_radiant_temperature_delta, 2) + \
                  1.94960053E-9 * air_temperature * np.power(air_velocity, 3) * np.power(mean_radiant_temperature_delta, 2) + \
                  -1.00361113E-8 * np.power(air_velocity, 4) * np.power(mean_radiant_temperature_delta, 2) + \
                  -1.21206673E-5 * np.power(mean_radiant_temperature_delta, 3) + \
                  -2.18203660E-7 * air_temperature * np.power(mean_radiant_temperature_delta, 3) + \
                  7.51269482E-9 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 3) + \
                  9.79063848E-11 * np.power(air_temperature, 3) * np.power(mean_radiant_temperature_delta, 3) + \
                  1.25006734E-6 * air_velocity * np.power(mean_radiant_temperature_delta, 3) + \
                  -1.81584736E-9 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 3) + \
                  -3.52197671E-10 * np.power(air_temperature, 2) * air_velocity * np.power(mean_radiant_temperature_delta, 3) + \
                  -3.36514630E-8 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 3) + \
                  1.35908359E-10 * air_temperature * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 3) + \
                  4.17032620E-10 * np.power(air_velocity, 3) * np.power(mean_radiant_temperature_delta, 3) + \
                  -1.30369025E-9 * np.power(mean_radiant_temperature_delta, 4) + \
                  4.13908461E-10 * air_temperature * np.power(mean_radiant_temperature_delta, 4) + \
                  9.22652254E-12 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 4) + \
                  -5.08220384E-9 * air_velocity * np.power(mean_radiant_temperature_delta, 4) + \
                  -2.24730961E-11 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 4) + \
                  1.17139133E-10 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 4) + \
                  6.62154879E-10 * np.power(mean_radiant_temperature_delta, 5) + \
                  4.03863260E-13 * air_temperature * np.power(mean_radiant_temperature_delta, 5) + \
                  1.95087203E-12 * air_velocity * np.power(mean_radiant_temperature_delta, 5) + \
                  -4.73602469E-12 * np.power(mean_radiant_temperature_delta, 6) + \
                  5.12733497 * saturated_vapor_pressure + \
                  -0.312788561 * air_temperature * saturated_vapor_pressure + \
                  -0.0196701861 * np.power(air_temperature, 2) * saturated_vapor_pressure + \
                  9.99690870E-4 * np.power(air_temperature, 3) * saturated_vapor_pressure + \
                  9.51738512E-6 * np.power(air_temperature, 4) * saturated_vapor_pressure + \
                  -4.66426341E-7 * np.power(air_temperature, 5) * saturated_vapor_pressure + \
                  0.548050612 * air_velocity * saturated_vapor_pressure + \
                  -0.00330552823 * air_temperature * air_velocity * saturated_vapor_pressure + \
                  -0.00164119440 * np.power(air_temperature, 2) * air_velocity * saturated_vapor_pressure + \
                  -5.16670694E-6 * np.power(air_temperature, 3) * air_velocity * saturated_vapor_pressure + \
                  9.52692432E-7 * np.power(air_temperature, 4) * air_velocity * saturated_vapor_pressure + \
                  -0.0429223622 * np.power(air_velocity, 2) * saturated_vapor_pressure + \
                  0.00500845667 * air_temperature * np.power(air_velocity, 2) * saturated_vapor_pressure + \
                  1.00601257E-6 * np.power(air_temperature, 2) * np.power(air_velocity, 2) * saturated_vapor_pressure + \
                  -1.81748644E-6 * np.power(air_temperature, 3) * np.power(air_velocity, 2) * saturated_vapor_pressure + \
                  -1.25813502E-3 * np.power(air_velocity, 3) * saturated_vapor_pressure + \
                  -1.79330391E-4 * air_temperature * np.power(air_velocity, 3) * saturated_vapor_pressure + \
                  2.34994441E-6 * np.power(air_temperature, 2) * np.power(air_velocity, 3) * saturated_vapor_pressure + \
                  1.29735808E-4 * np.power(air_velocity, 4) * saturated_vapor_pressure + \
                  1.29064870E-6 * air_temperature * np.power(air_velocity, 4) * saturated_vapor_pressure + \
                  -2.28558686E-6 * np.power(air_velocity, 5) * saturated_vapor_pressure + \
                  -0.0369476348 * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  0.00162325322 * air_temperature * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -3.14279680E-5 * np.power(air_temperature, 2) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  2.59835559E-6 * np.power(air_temperature, 3) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -4.77136523E-8 * np.power(air_temperature, 4) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  8.64203390E-3 * air_velocity * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -6.87405181E-4 * air_temperature * air_velocity * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -9.13863872E-6 * np.power(air_temperature, 2) * air_velocity * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  5.15916806E-7 * np.power(air_temperature, 3) * air_velocity * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -3.59217476E-5 * np.power(air_velocity, 2) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  3.28696511E-5 * air_temperature * np.power(air_velocity, 2) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -7.10542454E-7 * np.power(air_temperature, 2) * np.power(air_velocity, 2) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -1.24382300E-5 * np.power(air_velocity, 3) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -7.38584400E-9 * air_temperature * np.power(air_velocity, 3) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  2.20609296E-7 * np.power(air_velocity, 4) * mean_radiant_temperature_delta * saturated_vapor_pressure + \
                  -7.32469180E-4 * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -1.87381964E-5 * air_temperature * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  4.80925239E-6 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -8.75492040E-8 * np.power(air_temperature, 3) * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  2.77862930E-5 * air_velocity * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -5.06004592E-6 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  1.14325367E-7 * np.power(air_temperature, 2) * air_velocity * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  2.53016723E-6 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -1.72857035E-8 * air_temperature * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -3.95079398E-8 * np.power(air_velocity, 3) * np.power(mean_radiant_temperature_delta, 2) * saturated_vapor_pressure + \
                  -3.59413173E-7 * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  7.04388046E-7 * air_temperature * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  -1.89309167E-8 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  -4.79768731E-7 * air_velocity * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  7.96079978E-9 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  1.62897058E-9 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 3) * saturated_vapor_pressure + \
                  3.94367674E-8 * np.power(mean_radiant_temperature_delta, 4) * saturated_vapor_pressure + \
                  -1.18566247E-9 * air_temperature * np.power(mean_radiant_temperature_delta, 4) * saturated_vapor_pressure + \
                  3.34678041E-10 * air_velocity * np.power(mean_radiant_temperature_delta, 4) * saturated_vapor_pressure + \
                  -1.15606447E-10 * np.power(mean_radiant_temperature_delta, 5) * saturated_vapor_pressure + \
                  -2.80626406 * np.power(saturated_vapor_pressure, 2) + \
                  0.548712484 * air_temperature * np.power(saturated_vapor_pressure, 2) + \
                  -0.00399428410 * np.power(air_temperature, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -9.54009191E-4 * np.power(air_temperature, 3) * np.power(saturated_vapor_pressure, 2) + \
                  1.93090978E-5 * np.power(air_temperature, 4) * np.power(saturated_vapor_pressure, 2) + \
                  -0.308806365 * air_velocity * np.power(saturated_vapor_pressure, 2) + \
                  0.0116952364 * air_temperature * air_velocity * np.power(saturated_vapor_pressure, 2) + \
                  4.95271903E-4 * np.power(air_temperature, 2) * air_velocity * np.power(saturated_vapor_pressure, 2) + \
                  -1.90710882E-5 * np.power(air_temperature, 3) * air_velocity * np.power(saturated_vapor_pressure, 2) + \
                  0.00210787756 * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -6.98445738E-4 * air_temperature * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 2) + \
                  2.30109073E-5 * np.power(air_temperature, 2) * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 2) + \
                  4.17856590E-4 * np.power(air_velocity, 3) * np.power(saturated_vapor_pressure, 2) + \
                  -1.27043871E-5 * air_temperature * np.power(air_velocity, 3) * np.power(saturated_vapor_pressure, 2) + \
                  -3.04620472E-6 * np.power(air_velocity, 4) * np.power(saturated_vapor_pressure, 2) + \
                  0.0514507424 * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  -0.00432510997 * air_temperature * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  8.99281156E-5 * np.power(air_temperature, 2) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  -7.14663943E-7 * np.power(air_temperature, 3) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  -2.66016305E-4 * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  2.63789586E-4 * air_temperature * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  -7.01199003E-6 * np.power(air_temperature, 2) * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  -1.06823306E-4 * np.power(air_velocity, 2) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  3.61341136E-6 * air_temperature * np.power(air_velocity, 2) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  2.29748967E-7 * np.power(air_velocity, 3) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 2) + \
                  3.04788893E-4 * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -6.42070836E-5 * air_temperature * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  1.16257971E-6 * np.power(air_temperature, 2) * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  7.68023384E-6 * air_velocity * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -5.47446896E-7 * air_temperature * air_velocity * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -3.59937910E-8 * np.power(air_velocity, 2) * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 2) + \
                  -4.36497725E-6 * np.power(mean_radiant_temperature_delta, 3) * np.power(saturated_vapor_pressure, 2) + \
                  1.68737969E-7 * air_temperature * np.power(mean_radiant_temperature_delta, 3) * np.power(saturated_vapor_pressure, 2) + \
                  2.67489271E-8 * air_velocity * np.power(mean_radiant_temperature_delta, 3) * np.power(saturated_vapor_pressure, 2) + \
                  3.23926897E-9 * np.power(mean_radiant_temperature_delta, 4) * np.power(saturated_vapor_pressure, 2) + \
                  -0.0353874123 * np.power(saturated_vapor_pressure, 3) + \
                  -0.221201190 * air_temperature * np.power(saturated_vapor_pressure, 3) + \
                  0.0155126038 * np.power(air_temperature, 2) * np.power(saturated_vapor_pressure, 3) + \
                  -2.63917279E-4 * np.power(air_temperature, 3) * np.power(saturated_vapor_pressure, 3) + \
                  0.0453433455 * air_velocity * np.power(saturated_vapor_pressure, 3) + \
                  -0.00432943862 * air_temperature * air_velocity * np.power(saturated_vapor_pressure, 3) + \
                  1.45389826E-4 * np.power(air_temperature, 2) * air_velocity * np.power(saturated_vapor_pressure, 3) + \
                  2.17508610E-4 * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 3) + \
                  -6.66724702E-5 * air_temperature * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 3) + \
                  3.33217140E-5 * np.power(air_velocity, 3) * np.power(saturated_vapor_pressure, 3) + \
                  -0.00226921615 * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  3.80261982E-4 * air_temperature * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  -5.45314314E-9 * np.power(air_temperature, 2) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  -7.96355448E-4 * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  2.53458034E-5 * air_temperature * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  -6.31223658E-6 * np.power(air_velocity, 2) * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 3) + \
                  3.02122035E-4 * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 3) + \
                  -4.77403547E-6 * air_temperature * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 3) + \
                  1.73825715E-6 * air_velocity * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 3) + \
                  -4.09087898E-7 * np.power(mean_radiant_temperature_delta, 3) * np.power(saturated_vapor_pressure, 3) + \
                  0.614155345 * np.power(saturated_vapor_pressure, 4) + \
                  -0.0616755931 * air_temperature * np.power(saturated_vapor_pressure, 4) + \
                  0.00133374846 * np.power(air_temperature, 2) * np.power(saturated_vapor_pressure, 4) + \
                  0.00355375387 * air_velocity * np.power(saturated_vapor_pressure, 4) + \
                  -5.13027851E-4 * air_temperature * air_velocity * np.power(saturated_vapor_pressure, 4) + \
                  1.02449757E-4 * np.power(air_velocity, 2) * np.power(saturated_vapor_pressure, 4) + \
                  -0.00148526421 * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 4) + \
                  -4.11469183E-5 * air_temperature * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 4) + \
                  -6.80434415E-6 * air_velocity * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 4) + \
                  -9.77675906E-6 * np.power(mean_radiant_temperature_delta, 2) * np.power(saturated_vapor_pressure, 4) + \
                  0.0882773108 * np.power(saturated_vapor_pressure, 5) + \
                  -0.00301859306 * air_temperature * np.power(saturated_vapor_pressure, 5) + \
                  0.00104452989 * air_velocity * np.power(saturated_vapor_pressure, 5) + \
                  2.47090539E-4 * mean_radiant_temperature_delta * np.power(saturated_vapor_pressure, 5) + \
                  0.00148348065 * np.power(saturated_vapor_pressure, 6)

    return utci_approx[0]


def standard_effective_temperature(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity, metabolic_rate=1, clo_value=1):
    """
    Compute the standard effective temperature

    :param air_temperature: temperature (C)
    :type air_temperature
    :param mean_radiant_temperature: mean radiant temperature (C)
    :type mean_radiant_temperature
    :param air_velocity: wind speed (m/s)
    :type air_velocity
    :param relative_humidity: relative humidity (%)
    :type relative_humidity
    :param metabolic_rate: metabolic rate (met)
    :type metabolic_rate
    :param clo_value: clo value ()
    :type clo_value
    :return: standard effective temperature
    :rtype
    """
    # TODO: Add source documentation
    # TODO: Reformat docstring to NumPy style (https://docs.scipy.org/doc/numpy-1.15.0/docs/howto_document.html)
    # TODO: Add metabolic rates to docstring

    import math

    metabolic_rates = {
        "Sleeping": 0.7,
        "Reclining": 0.8,
        "Sitting": 1.0,
        "Typing": 1.1,
        "Standing": 1.2,
        "Driving": 1.5,
        "Cooking": 1.8,
        "Walking": 1.7,
        "Walking 2mph": 2.0,
        "Lifting 10lbs": 2.1,
        "Walking 3mph": 2.6,
        "House Cleaning": 2.7,
        "Basketball": 3,
        "Dancing": 3.4,
        "Walking 4mph": 3.8,
        "Lifting 100lbs": 4.0,
        "Shoveling": 4.4,
        "Running 9mph": 9.5
    }

    wme = 0  # External work done - usually around 0 - used for adding additional heat to person given uncontrollable factors

    # Function to find the saturation vapour pressure, used frequently throughout the comfPierceSET function.
    def findSaturatedVaporPressureTorr(T):
        # calculates Saturated Vapor Pressure (Torr) at Temperature T  (C)
        return np.exp(18.6686 - 4030.183 / (T + 235.0))

    # Key initial variables.
    VaporPressure = (relative_humidity * findSaturatedVaporPressureTorr(air_temperature)) / 100
    AirVelocity = max(air_velocity, 0.1)
    KCLO = 0.25
    BODYWEIGHT = 69.9
    BODYSURFACEAREA = 1.8258
    METFACTOR = 58.2
    SBC = 0.000000056697  # Stefan-Boltzmann constant (W/m2K4)
    CSW = 170
    CDIL = 120
    CSTR = 0.5

    TempSkinNeutral = 33.7  # setpoint (neutral) value for Tsk
    TempCoreNeutral = 36.8  # setpoint value for Tcr
    TempBodyNeutral = 36.49  # setpoint for Tb (.1*TempSkinNeutral + .9*TempCoreNeutral)
    SkinBloodFlowNeutral = 6.3  # neutral value for SkinBloodFlow

    # INITIAL VALUES - start of 1st experiment
    TempSkin = TempSkinNeutral
    TempCore = TempCoreNeutral
    SkinBloodFlow = SkinBloodFlowNeutral
    MSHIV = 0.0
    ALFA = 0.1
    ESK = 0.1 * metabolic_rate

    # Start new experiment here (for graded experiments)
    # UNIT CONVERSIONS (from input variables)

    p = 101325.0 / 1000  # This variable is the pressure of the atmosphere in kPa and was taken from the psychrometrics.js file of the CBE comfort tool.

    PressureInAtmospheres = p * 0.009869
    LTIME = 60
    TIMEH = LTIME / 60.0
    RCL = 0.155 * clo_value
    # AdjustICL(RCL, Conditions);  TH: I don't think this is used in the software

    FACL = 1.0 + 0.15 * clo_value  # % INCREASE IN BODY SURFACE AREA DUE TO CLOTHING
    LR = 2.2 / PressureInAtmospheres  # Lewis Relation is 2.2 at sea level
    RM = metabolic_rate * METFACTOR
    M = metabolic_rate * METFACTOR

    if clo_value <= 0:
        WCRIT = 0.38 * np.power(AirVelocity, -0.29)
        ICL = 1.0
    else:
        WCRIT = 0.59 * np.power(AirVelocity, -0.08)
        ICL = 0.45

    CHC = 3.0 * np.power(PressureInAtmospheres, 0.53)
    CHCV = 8.600001 * np.power((AirVelocity * PressureInAtmospheres), 0.53)
    CHC = max(CHC, CHCV)

    # initial estimate of Tcl
    CHR = 4.7
    CTC = CHR + CHC
    RA = 1.0 / (FACL * CTC)  # resistance of air layer to dry heat transfer
    TOP = (CHR * mean_radiant_temperature + CHC * air_temperature) / CTC
    TCL = TOP + (TempSkin - TOP) / (CTC * (RA + RCL))

    # ========================  BEGIN ITERATION
    #
    # Tcl and CHR are solved iteratively using: H(Tsk - To) = CTC(Tcl - To),
    # where H = 1/(Ra + Rcl) and Ra = 1/Facl*CTC

    TCL_OLD = TCL
    TIME = range(LTIME)
    flag = True
    for TIM in TIME:
        if flag == True:
            while abs(TCL - TCL_OLD) > 0.01:
                TCL_OLD = TCL
                CHR = 4.0 * SBC * np.power(((TCL + mean_radiant_temperature) / 2.0 + 273.15), 3.0) * 0.72
                CTC = CHR + CHC
                RA = 1.0 / (FACL * CTC)  # resistance of air layer to dry heat transfer
                TOP = (CHR * mean_radiant_temperature + CHC * air_temperature) / CTC
                TCL = (RA * TempSkin + RCL * TOP) / (RA + RCL)
        flag = False
        DRY = (TempSkin - TOP) / (RA + RCL)
        HFCS = (TempCore - TempSkin) * (5.28 + 1.163 * SkinBloodFlow)
        ERES = 0.0023 * M * (44.0 - VaporPressure)
        CRES = 0.0014 * M * (34.0 - air_temperature)
        SCR = M - HFCS - ERES - CRES - wme
        SSK = HFCS - DRY - ESK
        TCSK = 0.97 * ALFA * BODYWEIGHT
        TCCR = 0.97 * (1 - ALFA) * BODYWEIGHT
        DTSK = (SSK * BODYSURFACEAREA) / (TCSK * 60.0)  # //deg C per minute
        DTCR = SCR * BODYSURFACEAREA / (TCCR * 60.0)  # //deg C per minute
        TempSkin = TempSkin + DTSK
        TempCore = TempCore + DTCR
        TB = ALFA * TempSkin + (1 - ALFA) * TempCore
        SKSIG = TempSkin - TempSkinNeutral
        WARMS = (SKSIG > 0) * SKSIG
        COLDS = ((-1.0 * SKSIG) > 0) * (-1.0 * SKSIG)
        CRSIG = (TempCore - TempCoreNeutral)
        WARMC = (CRSIG > 0) * CRSIG
        COLDC = ((-1.0 * CRSIG) > 0) * (-1.0 * CRSIG)
        BDSIG = TB - TempBodyNeutral
        WARMB = (BDSIG > 0) * BDSIG
        COLDB = ((-1.0 * BDSIG) > 0) * (-1.0 * BDSIG)
        SkinBloodFlow = (SkinBloodFlowNeutral + CDIL * WARMC) / (1 + CSTR * COLDS)
        if SkinBloodFlow > 90.0: SkinBloodFlow = 90.0
        if SkinBloodFlow < 0.5: SkinBloodFlow = 0.5
        REGSW = CSW * WARMB * np.exp(WARMS / 10.7)
        if REGSW > 500.0: REGSW = 500.0
        ERSW = 0.68 * REGSW
        REA = 1.0 / (LR * FACL * CHC)  # evaporative resistance of air layer
        RECL = RCL / (LR * ICL)  # evaporative resistance of clothing (icl=.45)
        EMAX = (findSaturatedVaporPressureTorr(TempSkin) - VaporPressure) / (REA + RECL)
        PRSW = ERSW / EMAX
        PWET = 0.06 + 0.94 * PRSW
        EDIF = PWET * EMAX - ERSW
        ESK = ERSW + EDIF
        if PWET > WCRIT:
            PWET = WCRIT
            PRSW = WCRIT / 0.94
            ERSW = PRSW * EMAX
            EDIF = 0.06 * (1.0 - PRSW) * EMAX
            ESK = ERSW + EDIF
        if EMAX < 0:
            EDIF = 0
            ERSW = 0
            PWET = WCRIT
            PRSW = WCRIT
            ESK = EMAX
        ESK = ERSW + EDIF
        MSHIV = 19.4 * COLDS * COLDC
        M = RM + MSHIV
        ALFA = 0.0417737 + 0.7451833 / (SkinBloodFlow + .585417)

    # Define new heat flow terms, coeffs, and abbreviations
    STORE = M - wme - CRES - ERES - DRY - ESK  # rate of body heat storage
    HSK = DRY + ESK  # total heat loss from skin
    RN = M - wme  # net metabolic heat production
    ECOMF = 0.42 * (RN - (1 * METFACTOR))
    if ECOMF < 0.0: ECOMF = 0.0  # from Fanger
    EREQ = RN - ERES - CRES - DRY
    EMAX = EMAX * WCRIT
    HD = 1.0 / (RA + RCL)
    HE = 1.0 / (REA + RECL)
    W = PWET
    PSSK = findSaturatedVaporPressureTorr(TempSkin)
    # Definition of ASHRAE standard environment... denoted "S"
    CHRS = CHR
    if metabolic_rate < 0.85:
        CHCS = 3.0
    else:
        CHCS = 5.66 * np.power((metabolic_rate - 0.85), 0.39)
        if CHCS < 3.0: CHCS = 3.0

    CTCS = CHCS + CHRS
    RCLOS = 1.52 / ((metabolic_rate - wme / METFACTOR) + 0.6944) - 0.1835
    RCLS = 0.155 * RCLOS
    FACLS = 1.0 + KCLO * RCLOS
    FCLS = 1.0 / (1.0 + 0.155 * FACLS * CTCS * RCLOS)
    IMS = 0.45
    ICLS = IMS * CHCS / CTCS * (1 - FCLS) / (CHCS / CTCS - FCLS * IMS)
    RAS = 1.0 / (FACLS * CTCS)
    REAS = 1.0 / (LR * FACLS * CHCS)
    RECLS = RCLS / (LR * ICLS)
    HD_S = 1.0 / (RAS + RCLS)
    HE_S = 1.0 / (REAS + RECLS)

    # SET* (standardized humidity, clo_value, Pb, and CHC) determined using Newton's iterative solution
    # FNERRS is defined in the GENERAL SETUP section above

    DELTA = .0001
    dx = 100.0
    X_OLD = TempSkin - HSK / HD_S  # lower bound for SET
    while abs(dx) > .01:
        ERR1 = (HSK - HD_S * (TempSkin - X_OLD) - W * HE_S * (PSSK - 0.5 * findSaturatedVaporPressureTorr(X_OLD)))
        ERR2 = (HSK - HD_S * (TempSkin - (X_OLD + DELTA)) - W * HE_S * (
                PSSK - 0.5 * findSaturatedVaporPressureTorr((X_OLD + DELTA))))
        X = X_OLD - DELTA * ERR1 / (ERR2 - ERR1)
        dx = X - X_OLD
        X_OLD = X

    return X


def utci_openfield(self):
    universal_thermal_climate_index_openfield = universal_thermal_climate_index(
        self.dry_bulb_temperature,
        self.mean_radiant_temperature_openfield,
        self.pedestrian_wind_speed,
        self.relative_humidity
    )
    self.universal_thermal_climate_index_openfield = pd.Series(index=self.index, data=universal_thermal_climate_index_openfield, name="universal_thermal_climate_index_openfield")
    # self.df = pd.concat([self.df, self.universal_thermal_climate_index_openfield], axis=1).drop_duplicates(keep="last")
    print("Universal thermal climate index (openfield) calculations successful")
    return self.universal_thermal_climate_index_openfield


def utci_solar_adjusted(self):
    universal_thermal_climate_index_solar_adjusted = universal_thermal_climate_index(
        self.dry_bulb_temperature,
        self.mean_radiant_temperature_solar_adjusted,
        self.pedestrian_wind_speed,
        self.relative_humidity
    )
    self.universal_thermal_climate_index_solar_adjusted = pd.Series(index=self.index, data=universal_thermal_climate_index_solar_adjusted, name="universal_thermal_climate_index_solar_adjusted")
    # self.df = pd.concat([self.df, self.universal_thermal_climate_index_solar_adjusted], axis=1).drop_duplicates(keep="last")
    print("Universal thermal climate index (solar adjusted) calculations successful (SA)")
    return self.universal_thermal_climate_index_solar_adjusted
