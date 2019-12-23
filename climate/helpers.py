from climate.constants import *

import os
import math
import re
import numpy as np
from scipy import spatial


##########################
# Generic helper methods #
##########################


def closest_point(source_points, target_points, n_closest=1):
    """
    Find the closest n points from a set of source points to a set of target points
    :param source_points:
    :param target_points:
    :param n_closest:
    :return:
    """
    target_distances, target_indices = spatial.KDTree(target_points).query(source_points, n_closest)
    return target_distances, target_indices

def is_mac():
    if os.name == "posix":
        return True
    else:
        return False


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


def renamer(original):
    rename = {
        "dry_bulb_temperature": "Dry-Bulb Temperature",
        "dew_point_temperature": "Dew-Point Temperature",
        "relative_humidity": "Relative Humidity",
        "atmospheric_station_pressure": "Atmospheric Station Pressure",
        "extraterrestrial_horizontal_radiation": "Extraterrestrial Horizontal Radiation",
        "extraterrestrial_direct_normal_radiation": "Extraterrestrial Direct Normal Radiation",
        "horizontal_infrared_radiation_intensity": "Horizontal Infrared Radiation Intensity",
        "global_horizontal_radiation": "Global Horizontal Radiation",
        "direct_normal_radiation": "Direct Normal Radiation",
        "diffuse_horizontal_radiation": "Diffuse Horizontal Radiation",
        "global_horizontal_illuminance": "Global Horizontal Illuminance",
        "direct_normal_illuminance": "Direct Normal Illuminance",
        "diffuse_horizontal_illuminance": "Diffuse Horizontal Illuminance",
        "zenith_luminance": "Zenith Luminance",
        "wind_direction": "Wind Direction",
        "wind_speed": "Wind Speed",
        "total_sky_cover": "Total Sky Cover",
        "opaque_sky_cover": "Opaque Sky Cover",
        "visibility": "Visibility",
        "ceiling_height": "Ceiling Height",
        "present_weather_observation": "Present Weather Observation",
        "present_weather_codes": "Present Weather Codes",
        "precipitable_water": "Precipitable Water",
        "aerosol_optical_depth": "Aerosol Optical Depth",
        "snow_depth": "Snow Depth",
        "days_since_last_snowfall": "Days Since Last Snowfall",
        "albedo": "Albedo",
        "liquid_precipitation_depth": "Liquid Precipitation Depth",
        "liquid_precipitation_quantity": "Liquid Precipitation Quantity",
        "solar_apparent_zenith_angle": "Solar Apparent Zenith Angle",
        "solar_zenith_angle": "Solar Zenith Angle",
        "solar_apparent_elevation_angle": "Solar Apparent Elevation Angle",
        "solar_elevation_angle": "Solar Elevation Angle",
        "solar_azimuth_angle": "Solar Azimuth Angle",
        "solar_equation_of_time": "Solar Equation of Time",
        "humidity_ratio": "Humidity Ratio",
        "wet_bulb_temperature": "Wet Bulb Temperature",
        "partial_vapour_pressure_moist_air": "Partial Vapour Pressure of Moist air",
        "enthalpy": "Enthalpy",
        "specific_volume_moist_air": "Specific Volume of Moist Air",
        "degree_of_saturation": "Degree of Saturation",
    }
    return rename[original]


def radiation_rose_values(rose_vectors, sky_dome_patch_normal_vectors, sky_dome_patch_values):
    rose_vector_result = []
    for vec in rose_vectors:
        radiation = 0
        for patch_number, patch_vector in enumerate(sky_dome_patch_normal_vectors):
            vector_angle = angle_between(patch_vector, vec, degrees=True)
            if vector_angle < 90:
                radiation += sky_dome_patch_values[patch_number] * math.cos(math.radians(vector_angle))
        rose_vector_result.append(radiation)
    return np.array(rose_vector_result)


def wind_speed_at_height(source_wind_speed, source_wind_height, target_wind_height, terrain_roughness="Airport runway areas", log_method=True):
    """

    :param source_wind_speed: The wind speed to be translated (m/s)
    :param source_wind_height: The height at which the source wind speed was measured
    :param target_wind_height: The height to which the source wind speed will be translated
    :param terrain_roughness: A terrain roughness value from the European Wind Atlas, page 58, Figure 3.1: Roughness length, surface characteristics and roughness class.
    :param log_method:
    :return:
    """
    roughness = {
        "City": 1,
        "Forest": 0.8,
        "Suburbs": 0.5,
        "Shelter belts": 0.3,
        "Many trees and/or bushes": 0.2,
        "Farmland with closed appearance": 0.1,
        "Farmland with open appearance": 0.05,
        "Farmland with very few buildings, trees etc. airport areas with buildings and trees": 0.03,
        "Airport runway areas": 0.01,
        "Mown grass": 0.0075,
        "Bare soil (smooth)": 0.005,
        "Snow surfaces (smooth)": 0.001,
        "Sand surfaces (smooth)": 0.0003,
        "Water areas (lakes, fjords, open sea)": 0.0001
    }

    if source_wind_speed == 0:
        return 0
    if log_method:
        return source_wind_speed * (math.log(target_wind_height / roughness[terrain_roughness]) / math.log(source_wind_height / roughness[terrain_roughness]))
    wind_shear_exponent = 1 / 7
    return source_wind_speed * ((target_wind_height / source_wind_height) ** wind_shear_exponent)


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


def simple_combined_exterior_heat_transfer_coefficient(material_roughness, local_wind_speed):
    """
    The simple algorithm uses surface roughness and local surface windspeed to calculate the exterior heat transfer coefficient.

    Note that the simple correlation yields a combined convection and radiation heat transfer coefficient.

    Material roughness correlations are taken from Figure 1, Page 22.4, ASHRAE Handbook of Fundamentals (ASHRAE 1989).

    Source: http://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-020.html for further details.

    Parameters
    ----------
    material_roughness : string
        A surface roughness index. Acceptable values are: "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth" and "Very smooth"
    local_wind_speed : float
        Local wind speed calculated at the height above ground of the surface centroid

    Returns
    -------
    heat_transfer_coefficient
        The surface heat transfer coefficient for the variables passed
    """

    material = {
        "Very rough": {"D": 11.58, "E": 5.89, "F": 0},
        "Rough": {"D": 12.49, "E": 4.065, "F": 0.028},
        'Medium rough': {"D": 10.79, "E": 4.192, "F": 0},
        'Medium smooth': {"D": 8.23, "E": 4, "F": -0.057},
        'Smooth': {"D": 10.22, "E": 3.1, "F": 0},
        'Very smooth': {"D": 8.23, "E": 3.33, "F": -0.036},
    }

    exterior_heat_transfer_coefficient = material[material_roughness]["D"] + material[material_roughness]["E"] * local_wind_speed + material[material_roughness]["F"] * np.power(local_wind_speed, 2)
    return exterior_heat_transfer_coefficient
