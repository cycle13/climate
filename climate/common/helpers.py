# TODO - Add method to create mask for annual data given hours as list or hour as range, season or month

import warnings

from climate.common.constants import *

import pathlib
import os
import math
import re
import numpy as np
from scipy import spatial
from matplotlib.colors import LinearSegmentedColormap


def find_files(directory, extension=None):
    p = pathlib.Path(directory).glob('**/*')
    if extension is None:
        files = [x for x in p if x.is_file()]
    else:
        files = [x for x in p if x.is_file() and x.suffix == extension]
    return files


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


def chunk(enumerable, n=None, method="chunks"):
    """
    Partition a list into either n-chunks, or into a number of n-sized sub-lists

    Parameters
    ----------
    enumerable : ndarray
        A 1D array to split
    n : int
        Number of partitions or partition size to split list into
    method : str
        The method by which to split the array ("chunk" or "size")

    Returns
    -------
    chunked_list : list(lists)
        A list of sub-lists
    """

    chunked = None
    enumerable = np.array(enumerable)

    if method == "size":
        chunked = [list(enumerable[i:i + n]) for i in range(0, enumerable.shape[0], n)]

    elif method == "chunks":
        chunked = [list(i) for i in np.array_split(enumerable, n)]

    return chunked


def generate_directory(directory):
    """
    Create the directory passed

    Parameters
    ----------
    directory : str
        Directory to create

    Returns
    -------
    directory : str
        Created directory path
    """
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
    return directory


def hex_to_rgb(hex_string, normalise=False):
    """
    Convert an hexadecimal string to its RGB array equivalent

    Parameters
    ----------
    hex_string : str
        1D array containing RGB values in either 0-255 or 0-1 format

    Returns
    -------
    rgb_list : ndarray
        1D array containing RGB values in 0-255 format
    """
    hex_string = hex_string.lstrip('#')
    lv = len(hex_string)
    if normalise:
        return
    else:
        return np.array(tuple(int(hex_string[i:i + lv // 3], 16) for i in range(0, lv, lv // 3)))


def normalise_rgb(rgb_color):
    if sum([i < 1 for i in rgb_color]) > 1:
        warnings.warn("\n    Input color may be interpreted incorrectly as the composite RGB values are all less than 1. It's probably worth checking it's correct!", Warning)
        interpreted_color = np.interp(np.array(rgb_color), [0, 1], [0, 255]).tolist()
    else:
        interpreted_color = rgb_color
    return interpreted_color

def denormalise_rgb(rgb_color):
    if sum([i < 1 for i in rgb_color]) > 1:
        warnings.warn("\n    Input color may be interpreted incorrectly as the composite RGB values are all less than 1. It's probably worth checking it's correct!", Warning)
        interpreted_color = np.interp(np.array(rgb_color), [0, 255], [0, 1]).tolist()
    else:
        interpreted_color = rgb_color
    return interpreted_color


def rgb_to_hex(rgb_list):
    """
    Convert an RGB tuple/list to its hexadecimal equivalent

    Parameters
    ----------
    rgb_list : ndarray
        1D array containing RGB values in either 0-255 or 0-1 format

    Returns
    -------
    hex_color : str
        Hexadecimal representation of the input RGB color
    """
    a, b, c = normalise_rgb(rgb_list)
    return '#%02x%02x%02x' % (max(min([int(a), 255]), 0), max(min([int(b), 255]), 0), max(min([int(c), 255]), 0))


def gen_cmap(colors, n=100, r=False):
    translated_colors = []
    for color in colors:
        if type(color) == str:
            translated_colors.append(hex_to_rgb(color))
        else:
            translated_colors.append(normalise_rgb(color))
    cm = LinearSegmentedColormap.from_list('test', np.interp(translated_colors, [0, 255], [0, 1]), N=n)
    if r:
        return cm.reversed()
    else:
        return cm


####################
# Geometry methods #
####################


def unit_vector(vector):
    """
    Returns the angle between vectors 'v1' and 'v2'

    Parameters
    ----------
    vector : vector
        Non-unitized n-dimensional vector

    Returns
    -------
    unit_vector : ndarray
        The unitized form of the input vector
    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, degrees=False):
    """
    Returns the angle between vectors 'v1' and 'v2'

    Parameters
    ----------
    v1 : vector
        The first vector
    v2 : vector
        The second vector
    degrees : bool
        True for value returned in degrees, False for value returned in radians

    Returns
    -------
    angle : float
        The angle between the passed vectors
    """
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
        "dry_bulb_temperature": "Dry-Bulb Temperature (°C)",
        "dew_point_temperature": "Dew-Point Temperature (°C)",
        "relative_humidity": "Relative Humidity (%)",
        "atmospheric_station_pressure": "Atmospheric Station Pressure (Pa)",
        "extraterrestrial_horizontal_radiation": "Extraterrestrial Horizontal Radiation (W/m²)",
        "extraterrestrial_direct_normal_radiation": "Extraterrestrial Direct Normal Radiation (W/m²)",
        "horizontal_infrared_radiation_intensity": "Horizontal Infrared Radiation Intensity (W/m²)",
        "global_horizontal_radiation": "Global Horizontal Radiation (W/m²)",
        "direct_normal_radiation": "Direct Normal Radiation (W/m²)",
        "diffuse_horizontal_radiation": "Diffuse Horizontal Radiation (W/m²)",
        "global_horizontal_illuminance": "Global Horizontal Illuminance (lux)",
        "direct_normal_illuminance": "Direct Normal Illuminance (lux)",
        "diffuse_horizontal_illuminance": "Diffuse Horizontal Illuminance (lux)",
        "zenith_luminance": "Zenith Luminance (Cd/m²)",
        "wind_direction": "Wind Direction (degrees)",
        "wind_speed": "Wind Speed (m/s)",
        "total_sky_cover": "Total Sky Cover (tenths)",
        "opaque_sky_cover": "Opaque Sky Cover (tenths)",
        "visibility": "Visibility (km)",
        "ceiling_height": "Ceiling Height (m)",
        "present_weather_observation": "Present Weather Observation",
        "present_weather_codes": "Present Weather Codes",
        "precipitable_water": "Precipitable Water (mm)",
        "aerosol_optical_depth": "Aerosol Optical Depth (thousandths)",
        "snow_depth": "Snow Depth (cm)",
        "days_since_last_snowfall": "Days Since Last Snowfall",
        "albedo": "Albedo",
        "liquid_precipitation_depth": "Liquid Precipitation Depth (mm)",
        "liquid_precipitation_quantity": "Liquid Precipitation Quantity (hr)",
        "solar_apparent_zenith_angle": "Solar Apparent Zenith Angle (degrees)",
        "solar_zenith_angle": "Solar Zenith Angle (degrees)",
        "solar_apparent_elevation_angle": "Solar Apparent Elevation Angle (degrees)",
        "solar_elevation_angle": "Solar Elevation Angle (degrees)",
        "solar_azimuth_angle": "Solar Azimuth Angle (degrees)",
        "solar_equation_of_time": "Solar Equation of Time (minutes)",
        "humidity_ratio": "Humidity Ratio",
        "wet_bulb_temperature": "Wet Bulb Temperature (°C)",
        "partial_vapour_pressure_moist_air": "Partial Vapour Pressure of Moist air (Pa)",
        "enthalpy": "Enthalpy (J/kg)",
        "specific_volume_moist_air": "Specific Volume of Moist Air (m³/kg)",
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


def calc_sky_emissivity(dew_point_temperature, total_sky_cover):
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


def generate_numerous_vectors(samples=1000):
    """
    Create a set of rays from a point source

    Parameters
    ----------
    samples : int
        Number of rays to generate

    Returns
    -------
    vector : array(x, y, z)
        Vector direction for each ray
    theta : array(float)
        Angle (in radians) between point source plane and vector end-point
    sky : array(bool)
        denotes whether ray cast is above/below ground (True is above ground - towards the sky)
    """

    offset = 2 / samples
    increment = np.pi * (3 - np.sqrt(5))
    y = ((np.arange(samples) * offset) - 1) + (offset / 2)
    r = np.sqrt(1 - np.power(y, 2))
    phi = np.arange(samples) * increment
    x = np.cos(phi) * r
    z = np.sin(phi) * r
    vector = np.array([x, y, z]).T
    theta = np.fabs(np.arctan(z / np.sqrt(np.power(x, 2) + np.power(y, 2))))
    sky = z > 0
    return vector, theta, sky


def closest_point(source_points, target_points, n_closest=1):
    """
    Find the closest n-points within a set of target points from a set of source points.

    Parameters
    ----------
    source_points : array(x, y, z)
        Set of source points which will be used
    target_points : array(x, y, z)
        Set of target points which will be assessed for proximity
    n_closest : int
        The number of "near" points to return

    Returns
    -------
    target_distances : array(float)
        Distance between each source point and nearest target point/s
    target_indices : array(int)
        Indices of each target point near to the source point/s

    """
    target_distances, target_indices = spatial.KDTree(target_points).query(source_points, n_closest)
    return target_distances, target_indices


def sky_view_factor(altitude, is_sky):
    """
    Calculate the sky view factor between source and sample vector (including below horizon effects.

    Source: Jianxiang Huang, Jose Guillermo Cedeño-Laurent, John D. Spengler. 2014. CityComfort+: A simulation-based method for predicting mean radiant temperature in dense urban areas, https://doi.org/10.1016/j.buildenv.2014.05.019

    Parameters
    ----------
    altitude : float
        Angle between horizon and sample point (radians)
    is_sky : bool
        Different method used for below horizontal samples (accounting for person height)

    Returns
    -------
    sky_view_factor : array(float)
        Sky view factor value
    """
    az = np.pi / 4  # 45 degrees
    sky_view_factor = np.where(is_sky, 0.0355 * np.sin(altitude) + 2.33 * np.cos(altitude) * np.sqrt(
        0.0213 * np.power(np.cos(az), 2) + 0.0091 * np.power(np.sin(az), 2)), 0)
    return sky_view_factor


def radiation_on_square_metre(sample_vectors_altitude, sample_vectors_is_sky, sample_vector_radiation):
    A = sample_vectors_altitude[sample_vectors_is_sky]
    B = np.array([i[sample_vectors_is_sky] for i in sample_vector_radiation])
    C = np.sin(A) * B
    radiation = C.sum(axis=1)
    return radiation


def surface_temperature(dry_bulb_temperature, ground_temperature, ground_1m2_radiation, ground_emissivity,
                        ground_absorptivity, ground_k_value, ground_thickness,
                        ground_convective_heat_transfer_coefficient):
    a = ground_emissivity * STEFAN_BOLTZMANN_CONSTANT
    b = ground_k_value / ground_thickness + ground_convective_heat_transfer_coefficient
    c = -(ground_k_value * (
                ground_temperature + KELVIN) / ground_thickness + ground_1m2_radiation * ground_absorptivity + ground_convective_heat_transfer_coefficient * (
                      dry_bulb_temperature + KELVIN))

    Ts = []
    X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a)) + 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = -0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                               np.sqrt(3) * np.sqrt(
                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                               2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt((2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = 0.5 * np.sqrt((2 * b) / (a * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                     np.sqrt(3) * np.sqrt(
                                                 27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                     1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                              1 / 3)) - 1 / 2 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                                np.sqrt(3) * np.sqrt(
                                                            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                                1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))
    if X.sum() > 0:
        Ts.append(X)

    return Ts[0]