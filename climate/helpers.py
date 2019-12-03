import os
import math
import re
import numpy as np

##########################
# Generic helper methods #
##########################


def slugify(text):
    return re.sub(r'\W+', '', text)


def interpret(val):
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def chunk(enumerable, n_chunks):
    for i in range(0, len(enumerable), n_chunks):
        yield enumerable[i:i + n_chunks]


def generate_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def hex_to_rgb(hex_string):
    hex_string = hex_string.lstrip('#')
    lv = len(hex_string)
    return tuple(int(hex_string[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb_list):
    if (rgb_list[0] >= 1) | (rgb_list[1] >= 1) | (rgb_list[2] >= 1):
        return '#%02x%02x%02x' % (int(rgb_list[0]), int(rgb_list[1]), int(rgb_list[2]))
    else:
        return '#%02x%02x%02x' % (int(rgb_list[0] * 255), int(rgb_list[1] * 255), int(rgb_list[2] * 255))


####################
# Geometry methods #
####################

def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, degrees=False):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if degrees:
        return np.degrees(angle)
    else:
        return angle


###########################
# Generic climate methods #
###########################

def wind_speed_at_height(ws, h1, h2, rc=0, log=True):
    roughness = {
        0: 0.0002,  # Water surfaces: seas and Lakes
        0.2: 0.0005,  # Inlet water
        0.5: 0.0024,  # Open terrain with smooth surface, e.g. concrete, airport runways, mown grass etc.
        1: 0.03,  # Open agricultural land without fences and hedges; very gentle hills and the odd farmhouse
        1.5: 0.055,  # Agricultural land with a few buildings and 8 m high hedges separated by more than 1 km
        2: 0.1,  # Agricultural land with a few buildings and 8 m high hedges separated by approx. 500 m
        2.5: 0.2,  # Agricultural land with many trees, bushes and plants, or 8 m high hedges separated by approx. 250 m
        3: 0.4,  # Towns or villages with many or high hedges, forests and very rough and uneven terrain
        3.5: 0.6,  # Large towns with high buildings
        4: 1.6  # Large cities with high buildings and skyscrapers
    }
    if ws == 0:
        return 0
    if log:
        return ws * (math.log(h2 / roughness[rc]) / math.log(h1 / roughness[rc]))
    wind_shear_exponent = 1 / 7
    return ws * ((h2 / h1) ** wind_shear_exponent)


def ground_temperature_at_depth(depth, annual_average_temperature, annual_temperature_range, days_since_coldest_day, soil_diffusivity=0.01):
    soil_diffusivities = {
        "Rock": 0.02,
        "Wet clay": 0.015,
        "Wet sand": 0.01,
        "Dry clay": 0.002,
        "Dry sand": 0.001
    }

    w = 2 * math.pi / 365
    dd = math.sqrt(2 * soil_diffusivity / w)

    return annual_average_temperature - (annual_temperature_range / 2) * math.exp(-depth / dd) * math.cos((w * days_since_coldest_day) - (depth / dd))


#########################
# Psychrometric methods #
#########################


def celsius_to_kelvin(tc):
    return tc + 273.15


def partial_vapour_pressure(dry_bulb_temperature):
    c1 = -5674.5359
    c2 = 6.3925247
    c3 = -0.009677843
    c4 = 0.00000062215701
    c5 = 2.0747825E-09
    c6 = -9.484024E-13
    c7 = 4.1635019
    c8 = -5800.2206
    c9 = 1.3914993
    c10 = -0.048640239
    c11 = 0.000041764768
    c12 = -0.000000014452093
    c13 = 6.5459673

    t_k = celsius_to_kelvin(dry_bulb_temperature)  # convert Celsius into Kelvin

    if dry_bulb_temperature < 0:
        return math.exp(c1 / t_k + c2 + c3 * t_k + c4 * t_k ** 2 + c5 * t_k ** 3 + c6 * t_k ** 4 + c7 * math.log(t_k))
    else:
        return math.exp(c8 / t_k + c9 + c10 * t_k + c11 * t_k ** 2 + c12 * t_k ** 3 + c13 * math.log(t_k))


def humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
    # calculate partial water pressure
    pvp = relative_humidity * partial_vapour_pressure(dry_bulb_temperature)

    # calculate humidity ratio
    return 0.621945 * pvp / (pressure * 1000 - pvp)


def enthalpy_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
    return 1.005 * dry_bulb_temperature + humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure) * ( 2501 + 1.86 * dry_bulb_temperature)


def wet_bulb_temperature_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):

    # Initial Wet Bulb Temp (equal to dry_bulb_temperature)
    wet_bulb_temperature1 = dry_bulb_temperature

    # calculate saturation humidity ratio of wet bulb temperature
    saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))

    # calculate humidity ratio of partial vapour pressure
    humidity_ratio = humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure)

    # calculate wet bulb temperature
    wet_bulb_temperature2 = ((2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 4.186 * humidity_ratio - 2.326 * saturation_vapour_pressure)

    myerror = wet_bulb_temperature2 - wet_bulb_temperature1

    if dry_bulb_temperature < 0:
        while abs(myerror) > 0.001:
            wet_bulb_temperature1 = 0.1 * myerror + wet_bulb_temperature1
            saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
            wet_bulb_temperature2 = ((2830 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2830 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 2.1 * humidity_ratio - 0.24 * saturation_vapour_pressure)
            myerror = wet_bulb_temperature2 - wet_bulb_temperature1
    else:
        while abs(myerror) > 0.001:
            wet_bulb_temperature1 = 0.01 * myerror + wet_bulb_temperature1
            saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
            wet_bulb_temperature2 = ((2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 4.186 * humidity_ratio - 2.326 * saturation_vapour_pressure)
            myerror = wet_bulb_temperature2 - wet_bulb_temperature1

    return wet_bulb_temperature2


def temperature_enthalpy_humidity_ratio(enthalpy, humidity_ratio):
    """ Calculate temperature as a function of enthalpy and humidity ratio"""
    return (enthalpy - 2.5 * (humidity_ratio * 1000)) / (1.01 + (0.00189 * humidity_ratio * 1000))
