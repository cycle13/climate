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


###################
# NV method stuff #
###################



def thing1(Tin, Ta, emissivity, k, tickness, hc, Ein, absorptivity):
    d=5.67*(10**(-8))

    """
    Convective Heat Transfer for air
    http://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
    The convective heat transfer coefficient of air is approximately equal to
    hc = 10.45-v+10*v**(1/2)    (2)
    v = the relative speed of the object through the air (m/s)
    Note! - this is an empirical equation and can be used for velocities - v - from 2 to 20 m/s.
    """

    """
    Emissivity = Absorptivity at a particular wavelength (and direction for “non-diffuse” surfaces).
    In the case of a surface receiving radiation and emitting radiation there is no reason why these
    two values should be the same. This is because the incident radiation will be at one wavelength
    (or range of wavelengths), but the wavelength of emission depends on the temperature of the surface.
    """

    """
    The equation for convection can be expressed as:
    q = hc A dT
    q = heat transferred per unit time (W)
    A = heat transfer area of the surface (m2)
    hc= convective heat transfer coefficient of the process (W/(m2K) or W/(m2oC))
    dT = temperature difference between the surface and the bulk fluid (K or oC)

    Heat Transfer Coefficients - Units
        1 W/(m2K) = 0.85984 kcal/(h m2 oC) = 0.1761 Btu/(ft2 h oF)
        1 Btu/(ft2 h oF) = 5.678 W/(m2 K) = 4.882 kcal/(h m2 oC)
        1 kcal/(h m2 oC) = 1.163 W/(m2K) = 0.205 Btu/(ft2 h oF)

    http://www.engineersedge.com/heat_transfer/convection_heat_transfer.htm
    """

    Tin=Tin+273.15
    Ta=Ta+273.15
    a=emissivity*d
    b=k/tickness+hc
    c=-(k*Tin/tickness+Ein*absorptivity+hc*Ta)


    Ts=[]
    nnn=[]
    try:
        X = 1/2 * math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a))-1/2*math.sqrt(-(2*b)/(a *math.sqrt((4* (2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27 *a**2 *b**4-256* a**3* c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)))-(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)-(4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2* b**4-256* a**3*c**3)+9*a*b**2)**(1/3))
        if 0<X:
            Ts.append(X)
    except: nnn.append(1)

    try:
        X = 1/2 *math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a))+1/2*math.sqrt(-(2*b)/(a*math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)))-(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)-(4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3))
        if 0<X:
            Ts.append(X)
    except: nnn.append(2)

    try:
        X =-1/2 *math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a))-1/2*math.sqrt((2*b)/(a*math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)))-(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)-(4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3))
        if 0<X:
            Ts.append(X)
    except: nnn.append(3)

    try:
        X = 1/2 *math.sqrt((2*b)/(a*math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)))-(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a)-(4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3))-1/2*math.sqrt((4*(2/3)**(1/3)*c)/(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)+(math.sqrt(3)*math.sqrt(27*a**2*b**4-256*a**3*c**3)+9*a*b**2)**(1/3)/(2**(1/3)*3**(2/3)*a))
        if 0<X:
            Ts.append(X)
    except: nnn.append(4)
    print(nnn)

    Tsc=[]
    for i in Ts:
        Tsc.append(i-273)
    Eout=[]
    for i in Ts:
        E=a*(i**4)
        Eout.append(E)

def thing2(i, va):
    """
    Convective Heat Transfer for air
    http://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-020.html
    The roughness correlation is taken from Figure 1, Page 22.4, ASHRAE Handbook of Fundamentals (ASHRAE 1989).

    hc = D+E*va+F*va**2

    hc = heat transfer coefficient
    Va = local wind speed calculated at the height above ground of the surface centroid
    D, E, F = material roughness coefficients
    """
    Material = ["Stucco (Very Rough)", "Brick (Rough)", 'Concrete (Medium Rough)', 'Clear pine (Medium Smooth)',
                'Smooth Plaster(Smooth) ', 'Glass (Very Smooth)']
    D = [11.58, 12.49, 10.79, 8.23, 10.22, 8.23]
    E = [5.89, 4.065, 4.192, 4, 3.1, 3.33]
    F = [0, 0.028, 0, -0.057, 0, -0.036]

    hc = D[i] + E[i] * va + F[i] * va ** 2
    mt = Material[i]

def thing3(Ta, Skyemissivity):
    # http://bigladdersoftware.com/epx/docs/8-4/engineering-reference/climate-calculations.html
    σ = 5.6697 * 10 ** (-8)
    Ta = Ta + 273.15
    HorizontalIR = Skyemissivity * σ * (Ta ** 4)