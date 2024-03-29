# from datetime import datetime, timezone, timedelta
# import math
# import json
# import os
# import pandas as pd
# import numpy as np
# 
# def generate_directory(directory):
#     import os
#     if not os.path.exists(directory):
#         os.makedirs(directory)
#     return directory
# 
# 
# def interpret(val):
#     """Attempt to interpret the dtype for a single object
# 
#     Parameters
#     ----------
#     val : object
#         An unknown object of type int or float or str
# 
#     Returns
#     -------
#     out : object
#         The input object, reinterpreted as either an int, float or it's original dtype
# 
#     Examples
#     --------
#     >>> interpret("3.4")
#     3.4
#     >>> interpret("2")
#     2
#     >>> interpret("A1")
#     "A1"
# 
#     """
# 
#     try:
#         return int(val)
#     except ValueError:
#         try:
#             return float(val)
#         except ValueError:
#             return val
# 
# 
# bh_colors = {
#     "BH_LightRed": [188, 32, 75],
#     "BH_BrightRed": [213, 0, 50],
#     "BH_LightBlue": [141, 185, 202],
#     "BH_BrightGreen": [108, 194, 78],
#     "BH_Teal": [0, 164, 153],
#     "BH_Orange": [255, 143, 28],
#     "BH_DarkBlue": [0, 97, 127],
#     "BH_Green": [34, 136, 72],
#     "BH_Purple": [36, 19, 95],
#     "BH_DarkRed": [126, 45, 64],
#     "BH_Purple": [36, 19, 95],
#     "BH_Grey": [150, 140, 131],
#     "BH_BrightPurple": [112, 47, 138],
#     "BH_Pink": [231, 130, 169],
#     "BH_Purple": [36, 19, 95],
# }
# 
# def gen_cmap(colors, N=100, r=False):
#     from matplotlib.colors import LinearSegmentedColormap
#     cm = LinearSegmentedColormap.from_list('test', np.interp(colors, [0, 255], [0, 1]), N=N)
#     if r:
#         return cm.reversed()
#     else:
#         return cm
# 
# 
# 
# 
# def chunks(l, n):
#     """ Yield successive n-sized chunks from l
# 
#     Parameters
#     ----------
#     l : list
#         List of objects
#     n : int
#         Number of elements in each returned sub-list
# 
#     Returns
#     -------
#     generator
#         A generator object of l split into n-sized chunks
# 
#     Examples
#     --------
#     >>> a = [3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8]
#     >>> b = chunks(a, 3)
#     >>> b
#     <generator object chunks at 0x000001B078A8DE60>
#     >>> list(b)
#     [[3, 4, 5], [6, 7, 8], [9, 1, 2], [3, 4, 5], [6, 7, 8]]
#     >>> c = chunks(a, 6)
#     >>> list(c)
#     [[3, 4, 5, 6, 7, 8], [9, 1, 2, 3, 4, 5], [6, 7, 8]]
#     """
# 
#     for i in range(0, len(l), n):
#         yield l[i:i + n]
# 
# def celsius_to_kelvin(tc):
#     """
#     Convert temperature in Celsius to Kelvin
#     :param Tc: temperature in Celsius
#     :type Tc
#     :return: temperature in Kelvin
#     :rtype
#     """
#     return tc + 273.15
# 
# 
# def partial_vapor_pressure(ambient_pressure, humidity_ratio):
#     """ Function to compute partial vapor pressure in [kPa]
# 
#     From page 6.9 equation 38 in ASHRAE Fundamentals handbook (2005)
# 
#     :param ambient_pressure: kPa
#     :param humidity_ratio: kg/kg dry air
#     :return: Partial vapor pressure (kPa)
#     """
#     return ambient_pressure * humidity_ratio / (0.62198 + humidity_ratio)
# 
# 
# def saturation_vapor_pressure(dry_bulb_temperature):
#     ''' Function to compute saturation vapor pressure in [kPa]
# 
#     From ASHRAE Fundamentals handbook (2005) p 6.2, equations 5 and 6
#     Valid from -100C to 200 C
# 
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :return: Saturation vapor pressure (kPa)
#     '''
#     import math
# 
#     C1 = -5674.5359
#     C2 = 6.3925247
#     C3 = -0.009677843
#     C4 = 0.00000062215701
#     C5 = 2.0747825E-09
#     C6 = -9.484024E-13
#     C7 = 4.1635019
#     C8 = -5800.2206
#     C8 = -5800.2206
#     C9 = 1.3914993
#     C10 = -0.048640239
#     C11 = 0.000041764768
#     C12 = -0.000000014452093
#     C13 = 6.5459673
# 
#     TK = celsius_to_kelvin(dry_bulb_temperature)
# 
#     if TK <= 273.15:
#         result = math.exp(C1 / TK + C2 + C3 * TK + C4 * TK ** 2 + C5 * TK ** 3 +
#                           C6 * TK ** 4 + C7 * math.log(TK)) / 1000
#     else:
#         result = math.exp(C8 / TK + C9 + C10 * TK + C11 * TK ** 2 + C12 * TK ** 3 +
#                           C13 * math.log(TK)) / 1000
#     return result
# 
# 
# def humidity_ratio(dry_bulb_temperature, wet_bulb_temperature, ambient_pressure):
#     ''' Function to calculate humidity ratio [kg H2O/kg air]
# 
#     From ASHRAE Fundamentals handbook (2005)
# 
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :param wet_bulb_temperature: Wet-bulb temperature (C)
#     :param ambient_pressure: Ambient pressure (kPa)
#     :return: Humidity ratio (kg H2O/kg air)
#     '''
# 
#     Pws = saturation_vapor_pressure(wet_bulb_temperature)
#     Ws = 0.62198 * Pws / (ambient_pressure - Pws)  # Equation 23, p6.8
#     if dry_bulb_temperature >= 0:  # Equation 35, p6.9
#         result = (((2501 - 2.326 * wet_bulb_temperature) * Ws - 1.006 * (
#                 dry_bulb_temperature - wet_bulb_temperature)) / (
#                           2501 + 1.86 * dry_bulb_temperature - 4.186 * wet_bulb_temperature))
#     else:  # Equation 37, p6.9
#         result = (((2830 - 0.24 * wet_bulb_temperature) * Ws - 1.006 * (
#                 dry_bulb_temperature - wet_bulb_temperature)) / (
#                           2830 + 1.86 * dry_bulb_temperature - 2.1 * wet_bulb_temperature))
#     return result
# 
# 
# def relative_humidity1(dry_bulb_temperature, wet_bulb_temperature, ambient_pressure):
#     ''' Calculates relative humidity ratio [%]
# 
#     From ASHRAE Fundamentals handbook (2005)
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :param wet_bulb_temperature: Wet-bulb temperature (C)
#     :param ambient_pressure: Ambient pressure (kPa)
#     :return: Relative humidity (%)
# 
#     '''
# 
#     W = humidity_ratio(dry_bulb_temperature, wet_bulb_temperature, ambient_pressure)
#     result = partial_vapor_pressure(ambient_pressure, W) / saturation_vapor_pressure(
#         dry_bulb_temperature)  # Equation 24, p6.8
#     return result
# 
# 
# def relative_humidity2(dry_bulb_temperature, ambient_pressure, humidity_ratio):
#     ''' Calculates the relative humidity ratio [%]
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :param ambient_pressure: Ambient pressure (kPa)
#     :param humidity_ratio: Humidity ratio (kg H2O/kg air)
#     :return: Relative humidity (%)
# 
#     '''
# 
#     Pw = partial_vapor_pressure(ambient_pressure, humidity_ratio)
#     Pws = saturation_vapor_pressure(dry_bulb_temperature)
#     result = Pw / Pws
#     return result
# 
# 
# def wet_bulb_temperature(dry_bulb_temperature, relative_humidity, ambient_pressure):
#     ''' Calculates the Wet Bulb temperature [C]
# 
#     TODO: Add source in ASHARE Handbook
#     Uses Newton-Rhapson iteration to converge quickly
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param ambient_pressure: Ambient pressure (kPa)
#     :return: Wet-bulb temperature (C)
# 
#     '''
# 
#     W_normal = relative_humidity2(dry_bulb_temperature, relative_humidity, ambient_pressure)
#     result = dry_bulb_temperature
# 
#     # Solves to within 0.001% accuracy using Newton-Rhapson
#     W_new = relative_humidity1(dry_bulb_temperature, result, ambient_pressure)
#     while abs((W_new - W_normal) / W_normal) > 0.00001:
#         W_new2 = relative_humidity1(dry_bulb_temperature, result - 0.001, ambient_pressure)
#         dw_dtwb = (W_new - W_new2) / 0.001
#         result = result - (W_new - W_normal) / dw_dtwb
#         W_new = relative_humidity1(dry_bulb_temperature, result, ambient_pressure)
#     return result
# 
# 
# def specific_enthalpy(dry_bulb_temperature, humidity_ratio):
#     ''' Calculates specific_enthalpy in kJ/kg (dry air)
# 
#         From from 2005 ASHRAE Handbook - Fundamentals - SI P6.9 eqn 32
# 
#         :param dry_bulb_temperature: Dry-bulb temperature (C)
#         :param humidity_ratio: Humidity ratio (kg H2O/kg air)
#         :return: specific_enthalpy (kJ/kg dry air)
# 
#     '''
# 
#     result = 1.006 * dry_bulb_temperature + humidity_ratio * (2501 + 1.86 * dry_bulb_temperature)
#     return result
# 
# 
# def dry_bulb_temperature(specific_enthalpy, humidity_ratio):
#     ''' Calculates Dry Bulb Temperature [C]
# 
#     Warning 0 state is 0C, 0%RH and 1ATM
# 
#     :param humidity_ratio: Humidity ratio (kg H2O/kg air)
#     :param: specific_enthalpy (kJ/kg dry air)
#     :return: Dry-bulb temperature (C)
# 
#     '''
#     result = (specific_enthalpy - (2501 * humidity_ratio)) / (1.006 + (1.86 * humidity_ratio))
#     return result
# 
# 
# def dew_point_temperature2(ambient_pressure, humidity_ratio):
#     ''' Function to compute the dew point temperature (deg C)
# 
#     From page 6.9 equation 39 and 40 in ASHRAE Fundamentals handbook (2005)
#     Valid for Dew Points less than 93C
# 
#     :param ambient_pressure: Ambient pressure (kPa)
#     :param humidity_ratio: Humidity ratio (kg H2O/kg air)
#     :return: Dew-point temperature (C)
# 
#     '''
# 
#     import math
# 
#     C14 = 6.54
#     C15 = 14.526
#     C16 = 0.7389
#     C17 = 0.09486
#     C18 = 0.4569
# 
#     Pw = partial_vapor_pressure(ambient_pressure, humidity_ratio)
#     alpha = math.log(Pw)
#     Tdp1 = C14 + C15 * alpha + C16 * alpha ** 2 + C17 * alpha ** 3 + C18 * Pw ** 0.1984
#     Tdp2 = 6.09 + 12.608 * alpha + 0.4959 * alpha ** 2
#     if Tdp1 >= 0:
#         result = Tdp1
#     else:
#         result = Tdp2
#     return result
# 
# 
# def dry_air_density(ambient_pressure, dry_bulb_temperature, humidity_ratio):
#     ''' Function to compute the dry air density [kg/m**3]
# 
#     From page 6.8 equation 28 ASHRAE Fundamentals handbook (2005)
#     Note that total density of air-h2o mixture is:
#     rho_air_h2o = rho_dry_air * (1 + humidity_ratio)
#     gas constant for dry air
# 
#     :param ambient_pressure: Ambient pressure (kPa)
#     :param dry_bulb_temperature: Dry-bulb temperature (C)
#     :param humidity_ratio: Humidity ratio (kg H2O/kg air)
#     :return: Dry air density (kg/m**3)
#     '''
# 
#     rho_dry_air = 287.055
#     result = 1000 * ambient_pressure / (
#                 rho_dry_air * (celsius_to_kelvin(dry_bulb_temperature)) * (1 + 1.6078 * humidity_ratio))
#     return result
# 
# 
# def solar_altitude(latitude, longitude, year, month, day, hour, time_zone, minute=0, second=0):
#     """
#     Determine suns altitude in sky from location on ground and timezone aware datetime
#     :param latitude: Latitude, between -90 and +90 (distance from equator, +Ve is north and -Ve is south)
#     :type latitude
#     :param longitude: Longitude, between -180 and +180 (distance from prime meridian, +Ve is east and -Ve is west)
#     :type longitude
#     :param year: Year
#     :type year: int
#     :param month: Month
#     :type month: int
#     :param day: Day
#     :type day: int
#     :param hour: Hour
#     :type hour: int
#     :param minute: Minute
#     :type minute: int
#     :param second: Second
#     :type second: int
#     :param time_zone: Number of hours offset from UTC between -12 and +12
#     :type time_zone:
#     :return: solar altitude
#     :rtype: (float)
#     # TODO: Replace with standard python method
#     """
#     from pysolar.solar import get_altitude_fast
#     from datetime import datetime, timezone, timedelta
# 
#     date_time = datetime(year, month, day, hour, minute, second, 0)
#     date_time = date_time.replace(tzinfo=timezone(timedelta(hours=time_zone)))
# 
#     altitude = get_altitude_fast(latitude, longitude, date_time)
# 
#     return altitude
# 
# 
# def solar_azimuth(latitude, longitude, year, month, day, hour, time_zone, minute=0, second=0):
#     """
#     Determine suns azimuth in sky from location on ground and timezone aware datetime
#     :param latitude: Latitude, between -90 and +90 (distance from equator, +Ve is north and -Ve is south)
#     :type latitude
#     :param longitude: Longitude, between -180 and +180 (distance from prime meridian, +Ve is east and -Ve is west)
#     :type longitude
#     :param year: Year
#     :type year: int
#     :param month: Month
#     :type month: int
#     :param day: Day
#     :type day: int
#     :param hour: Hour
#     :type hour: int
#     :param minute: Minute
#     :type minute: int
#     :param second: Second
#     :type second: int
#     :param time_zone: Number of hours offset from UTC between -12 and +12
#     :type time_zone:
#     :return: solar altitude
#     :rtype: (float)
#     # TODO: Replace with standard python method
#     """
#     from pysolar.solar import get_azimuth_fast
#     from datetime import datetime, timezone, timedelta
# 
#     date_time = datetime(year, month, day, hour, minute, second, 0)
#     date_time = date_time.replace(tzinfo=timezone(timedelta(hours=time_zone)))
# 
#     azimuth = get_azimuth_fast(latitude, longitude, date_time)
# 
#     return azimuth
# 
# 
# def convert_altitude(altitude):
#     """
#     Converts an altitude (from 0 below, to 180 above), to a range between 0 at horizon to 90 at directly above
#     :param altitude: Original altitude value
#     :return altitude_adjusted: Adjusted altitude value
#     """
#     if altitude < 0:
#         altitude_adjusted = 0
#     elif altitude > 90:
#         altitude_adjusted = 90 - altitude
#     else:
#         altitude_adjusted = altitude
#     return altitude_adjusted
# 
# 
# def convert_azimuth(azimuth):
#     """
#     Converts an azimuth (from 0 at North, moving clockwise), to a range between 0 at North and 180 South (both clock and anti-clockwise)
#     :param azimuth: Original azimuth value
#     :return azimuth_adjusted: Adjusted azimuth value
#     """
#     if azimuth < 180:
#         azimuth_adjusted = azimuth
#     elif azimuth > 180:
#         azimuth_adjusted = 360 - azimuth
#     return azimuth_adjusted
# 
# 
# def fanger_solar_exposure(azimuth, altitude, posture=1):
#     """
#     Determine the exposure factor based on human posture and sun position
#     :param azimuth: Solar azimuth, between 0 and 180 (inclusive)
#     :type azimuth
#     :param altitude: Solar altitude, between 0 and 90 (inclusive)
#     :type altitude
#     :param posture: Human posture. 0=Sitting, 1=Standing
#     :type posture: int
#     :return: solar exposure
#     :rtype
#     """
#     from scipy.interpolate import bisplev
# 
#     # Following lines commented dealt with sun position range issues. convert_azimuth() and convert()_altitude now handle this.
#     # if azimuth < 0:
#     #     raise ValueError('Azimuth is outside the expected range of 0 to 180.')
#     # elif azimuth > 180:
#     #     raise ValueError('Azimuth is outside the expected range of 0 to 180.')
#     # elif altitude < 0:
#     #     raise ValueError('Altitude is outside the expected range of 0 to 90.')
#     # elif altitude > 90:
#     #     raise ValueError('Altitude is outside the expected range of 0 to 90.')
# 
#     # Sitting bivariate B-spline and its derivatives
#     sitting_bisplrep = [[0, 0, 0, 0, 0, 180, 180, 180, 180, 180], [0, 0, 0, 0, 0, 90, 90, 90, 90, 90],
#                         [0.42999558258469306, 0.49777802226651985, 0.2858541264803382, 0.2636331839991635,
#                          0.10059901058304405, 0.5904998653021177, 0.6393605287969937, 0.41177803047742195,
#                          0.16397939147762605, 0.1145630272512949, -0.07290688451711066, -0.0877360565501316,
#                          0.03719879969518147, 0.06771059788029093, 0.09444998526069391, 0.5573351684449549,
#                          0.6212235986152396, 0.3384990152297299, 0.25505266892999545, 0.1011441730110879,
#                          0.4432956571996788, 0.49809124858382825, 0.29471168936411446, 0.19682482035937438,
#                          0.10008130856803796], 4, 4]
# 
#     # Standing bivariate B-spline and its derivatives
#     standing_bisplrep = [[0, 0, 0, 0, 0, 180, 180, 180, 180, 180], [0, 0, 0, 0, 0, 90, 90, 90, 90, 90],
#                          [0.365433469329803, 0.41471995039390336, 0.3539202584010255, 0.35205668670475776,
#                           0.21505967838173534, 0.5304745700779437, 0.6180584137541132, 0.11434859278048302,
#                           0.4862162611010728, 0.20252438358996272, -0.015147290187610778, 0.22189948439503024,
#                           0.6990946114268216, -0.000718703369787728, 0.22472889635480628, 0.5176922764465676,
#                           0.35055123160310636, -0.0032935618498728487, 0.3404006983313149, 0.19936403473400507,
#                           0.37870178660536147, 0.24613731159172733, 0.06300314787643235, 0.23364607863218287,
#                           0.2171651821703637], 4, 4]
# 
#     if posture == 0:
#         solar_exposure = bisplev(convert_azimuth(azimuth), convert_altitude(altitude), sitting_bisplrep)
#     elif posture == 1:
#         solar_exposure = bisplev(convert_azimuth(azimuth), convert_altitude(altitude), standing_bisplrep)
# 
#     return solar_exposure
# 
# 
# def solar_adjusted_mean_radiant_temperature(dry_bulb_temperature, diffuse_horizontal_radiation,
#                                             global_horizontal_radiation, direct_normal_radiation,
#                                             solar_azimuth, solar_altitude,
#                                             clothing_absorbtivity = 0.7, ground_reflectivity = 0.4,
#                                             shading_transmissivity = 1.0):
#     """
#     This function is based on the ladybug method by Mostapha Roudsari, building on formulas translating solar radiation
#     into an effective radiant field and solar-adjusted mean radiant temperature. For further details on this
#     approximation, see Arens, Edward; Huang, Li; Hoyt, Tyler; Zhou, Xin; Shiavon, Stefano. (2014). Modeling the comfort
#     effects of short-wave solar radiation indoors.  Indoor Environmental Quality (IEQ). http://escholarship.org/uc/item/89m1h2dg#page-4
#     :param dry_bulb_temperature: Dry-bulb temperature
#     :param diffuse_horizontal_radiation: Diffuse solar horizontal radiation
#     :param global_horizontal_radiation: Global solar horizontal radiation
#     :param direct_normal_radiation: Direct solar normal radiation
#     :param solar_azimuth: Solar azimuth
#     :param solar_altitude: Solar altitude
#     :param clothing_absorbtivity: The fraction of solar radiation absorbed by the human body. The default is set to 0.7 for (average/brown) skin and average clothing. Increase this value for darker skin or darker clothing.
#     :param ground_reflectivity:  The fraction of solar radiation reflected off of the ground. By default, this is set to 0.4, characteristic of concrete paving. 0.25 would correspond to outdoor grass or dry bare soil.
#     :param shading_transmissivity: The fraction of shading transmissivity. exposed=1, shaded=0.
#     :return:
#     """
#     # TODO: Add error handling, range and type checking
#     if solar_altitude <= 0:
#         solar_adjusted_mrt = dry_bulb_temperature
#     else:
#         fraction_body_exposed_radiation = 0.71  # Fraction of the body visible to radiation. 0.725 is for standing, 0.696 for sitting, 0.68 for lying down. A nominal value of 0.71 is used.
#         radiative_heat_transfer_coefficient = 6.012  # A good guess at the radiative heat transfer coefficient
#         projected_area_fraction = fanger_solar_exposure(solar_azimuth, solar_altitude, posture=1)
#         effective_radiant_flux = ((0.5 * fraction_body_exposed_radiation * (
#                 diffuse_horizontal_radiation + (global_horizontal_radiation * ground_reflectivity)) + (
#                                            fraction_body_exposed_radiation * projected_area_fraction * direct_normal_radiation)) * shading_transmissivity) * (
#                                          clothing_absorbtivity / 0.95)
#         mean_radiant_temperature_delta = (
#                 effective_radiant_flux / (fraction_body_exposed_radiation * radiative_heat_transfer_coefficient))
#         solar_adjusted_mrt = dry_bulb_temperature + mean_radiant_temperature_delta
#     return solar_adjusted_mrt
# 
# 
# def wind_speed_at_height(ws, h1, h2, rc=0, log=True):
#     """
#     Convert wind speed from one height to another using the specific method
#     :param ws: wind speed (m/s)
#     :type ws
#     :param h1: original height (m)
#     :type h1
#     :param h2: target height (m)
#     :type h2
#     :param rc: roughness class corresponding to length of obstacles
#     :type rc
#     :param log: method of calculation - either log or power
#     :type log: bool
#     :return: temperature in Fahrenheit
#     :rtype
#     """
# 
#     import math
# 
#     roughness = {
#         0: 0.0002,  # Water surfaces: seas and Lakes
#         0.2: 0.0005,  # Inlet water
#         0.5: 0.0024,  # Open terrain with smooth surface, e.g. concrete, airport runways, mown grass etc.
#         1: 0.03,  # Open agricultural land without fences and hedges; very gentle hills and the odd farmhouse
#         1.5: 0.055,  # Agricultural land with a few buildings and 8 m high hedges separated by more than 1 km
#         2: 0.1,  # Agricultural land with a few buildings and 8 m high hedges separated by approx. 500 m
#         2.5: 0.2,  # Agricultural land with many trees, bushes and plants, or 8 m high hedges separated by approx. 250 m
#         3: 0.4,  # Towns or villages with many or high hedges, forests and very rough and uneven terrain
#         3.5: 0.6,  # Large towns with high buildings
#         4: 1.6  # Large cities with high buildings and skyscrapers
#     }
#     if ws == 0:
#         return 0
#     if log:
#         return ws * (math.log(h2 / roughness[rc]) / math.log(h1 / roughness[rc]))
#     wind_shear_exponent = 1 / 7
#     return ws * ((h2 / h1) ** wind_shear_exponent)
# 
# 
# def ground_temperature(depth, annual_average_temperature, annual_temperature_range, days_since_coldest_day,
#                        soil_diffusivity=0.01):
#     """
#     Estimate ground temperature at a depth from surfaces based on air temperature
#     :param depth: depth (m)
#     :type depth
#     :param annual_average_temperature: average annual air temperature (C)
#     :type annual_average_temperature
#     :param annual_temperature_range: difference between min and max annual air temperatures (C)
#     :type annual_temperature_range
#     :param days_since_coldest_day: days since the coldest day in the year
#     :type days_since_coldest_day: int
#     :param soil_diffusivity: soil diffusivity (cm2/s)
#     :type soil_diffusivity
# 
#     soil_diffusivities = {
#         "Rock": 0.02,
#         "Wet clay": 0.015,
#         "Wet sand": 0.01,
#         "Dry clay": 0.002,
#         "Dry sand": 0.001
#     }
#     """
#     import math
# 
#     w = 2 * math.pi / 365
#     dd = math.sqrt(2 * soil_diffusivity / w)
# 
#     return annual_average_temperature - (annual_temperature_range / 2) * math.exp(-depth / dd) * math.cos(
#         (w * days_since_coldest_day) - (depth / dd))
# 
# 
# def celsius_to_fahrenheit(Tc):
#     """
#     Convert temperature in Celsius to Fahrenheit
#     :param Tc: temperature in Celsius
#     :type Tc
#     :return: temperature in Fahrenheit
#     :rtype
#     """
#     Tf = (Tc * (9 / 5)) + 32
#     return Tf
# 
# 
# def fahrenheit_to_celsius(Tf):
#     """
#     Convert temperature in Fahrenheit to Celsius
#     :param Tf: temperature in Fahrenheit
#     :type Tf
#     :return: temperature in Celsius
#     :rtype
#     """
#     Tc = (Tf - 32) * (5 / 9)
#     return Tc
# 
# 
# def celsius_to_kelvin(Tc):
#     """
#     Convert temperature in Celsius to Kelvin
#     :param Tc: temperature in Celsius
#     :type Tc
#     :return: temperature in Kelvin
#     :rtype
#     """
#     Tk = Tc + 273.15
#     return Tk
# 
# 
# def dewpoint_temperature1(Ta, rh):
#     # TODO: Add check that rh is between 0.0 and 100.0 for validity - if not this may not be as accurate!
#     """
#     Calculate dew-point temperature from air temperature and relative humidity using MET4 and MET4A method
#     for Tc in range:  0C < Tc < 60C and Tdp in range:  0C < Tdp < 50C, uncertainty is +/- 0.4C
#     :param Ta: temperature (Celsius)
#     :type Ta
#     :param rh: relative humidity (%)
#     :type rh
#     :return: dew-point temperature (Celsius)
#     :rtype
#     """
# 
#     import math
#     import logging
# 
#     if not 0 <= rh <= 100:
#         logging.warning(" Value input for RH is outside the range suitable for accurate Dewpoint estimation\n")
# 
#     a = 17.27
#     b = 237.7
#     rh = rh / 100
#     Tdp = (b * (((a * Ta) / (b + Ta)) + math.log(rh))) / (a - (((a * Ta) / (b + Ta)) + math.log(rh)))
# 
#     return Tdp
# 
# 
# def humidity_index(Ta, Tdp):
#     """
#     Calulate the humidity index (humidex) from air temperature and relative humidity -
#     https://web.archive.org/web/20130627223738/http://climate.weatheroffice.gc.ca/prods_servs/normals_documentation_e.html
#     :param Ta: temperature (Celcius)
#     :type Ta
#     :param Tdp: dew-point temperature (Celsius)
#     :type Tdp
#     :return: humidity index (%), humidity index category (0-4), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     import math
# 
#     dewpoint_k = celsius_to_kelvin(Tdp)
# 
#     e = 6.11 * math.exp(5417.7530 * ((1 / 273.16) - (1 / dewpoint_k)))
#     h = (0.5555) * (e - 10.0)
#     humidity_index = Ta + h
# 
#     humidity_index_categories = {
#         0: "Little to no discomfort",
#         1: "Some discomfort",
#         2: "Great discomfort; avoid exertion",
#         3: "Dangerous; heat stroke quite possible"
#     }
# 
#     humidity_index_category = None
#     comfortable = None
# 
#     if humidity_index < 30:
#         # Little to no discomfort
#         humidity_index_category = 0
#         comfortable = True
#     elif 30 <= humidity_index < 39:
#         # Some discomfort
#         humidity_index_category = 1
#         comfortable = False
#     elif 40 <= humidity_index < 45:
#         # Great discomfort; avoid exertion
#         humidity_index_category = 2
#         comfortable = False
#     elif humidity_index >= 45:
#         # Dangerous; heat stroke quite possible
#         humidity_index_category = 3
#         comfortable = False
# 
#     return humidity_index, humidity_index_category, comfortable
# 
# 
# def heat_index(Ta, rh):
#     """
#     Compute the Heat Index in accordance with the regression defined in the National Weather Service (NWS) Technical
#     Attachment (SR 90-23) - for details see https://www.weather.gov/safety/heat-index
#     :param Ta: temperature in Celcius
#     :type Ta
#     :param rh: relative humidity
#     :type rh
#     :return: heat index in Celcius, heat index category (0-4), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     import math
# 
#     Tf = celsius_to_fahrenheit(Ta)
# 
#     if Tf < 80:
#         heat_index_fahrenheit = 0.5 * (Tf + 61.0 + ((Tf - 68.0) * 1.2) + (rh * 0.094))
#     else:
#         heat_index_fahrenheit = -42.379 + 2.04901523 * Tf + 10.14333127 * rh - 0.22475541 * Tf * rh - 6.83783 * (
#                 10 ** (-3)) * (Tf ** (2)) - 5.481717 * (10 ** (-2)) * (rh ** (2)) + 1.22874 * (10 ** (-3)) * (
#                                         Tf ** (2)) * (rh) + 8.5282 * (10 ** (-4)) * (Tf) * (rh ** (2)) - 1.99 * (
#                                         10 ** (-6)) * (Tf ** (2)) * (rh ** (2))
#         if (Tf >= 80) and (Tf <= 112) and (rh < 13):
#             adjust = ((13 - rh) / 4) * math.sqrt((17 - abs(Tf - 95)) / 17)
#             heat_index_fahrenheit = heat_index_fahrenheit - adjust
#         elif (Tf >= 80) and (Tf <= 87) and (rh > 85):
#             adjust = ((rh - 85) / 10) * ((87 - Tf) / 5)
#             heat_index_fahrenheit = heat_index_fahrenheit + adjust
# 
#     heat_index_categories = {
#         0: None,
#         1: "Caution",
#         2: "Extreme caution",
#         3: "Danger",
#         4: "Extreme danger"
#     }
#     heat_index_category = None
#     comfortable = None
# 
#     if heat_index_fahrenheit < 80:
#         heat_index_category = 0
#         comfortable = True
#     elif 80 <= heat_index_fahrenheit < 90:
#         heat_index_category = 1
#         comfortable = False
#     elif 90 <= heat_index_fahrenheit < 105:
#         heat_index_category = 2
#         comfortable = False
#     elif 105 <= heat_index_fahrenheit < 130:
#         heat_index_category = 3
#         comfortable = False
#     elif heat_index_fahrenheit >= 130:
#         heat_index_category = 4
#         comfortable = False
# 
#     heat_index = fahrenheit_to_celsius(heat_index_fahrenheit)
# 
#     return heat_index, heat_index_category, comfortable
# 
# 
# def discomfort_index(Ta, rh):
#     """
#     Compute the Heat Index (Thom's Index) in accordance with the regression defined in Thom, E.C. (1959):
#     The discomfort index. Weather wise, 12: 5760
#     :param Ta: temperature in Celcius
#     :type Ta
#     :param rh: relative humidity
#     :type rh
#     :return: discomfort index in Celcius, discomfort index category (0-9), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     discomfort_index = Ta - (0.55 - 0.0055 * rh) * (Ta - 14.5)
# 
#     discomfort_index_categories = {
#         0: "Hyperglacial",
#         1: "Glacial",
#         2: "Extremely cold",
#         3: "Very cold",
#         4: "Cold",
#         5: "Cool",
#         6: "Comfortable",
#         7: "Hot",
#         8: "Very hot",
#         9: "Torrid"
#     }
# 
#     # Categories from Kyle WJ (1994) The human bioclimate of Hong-Kong.
#     discomfort_index_category = None
#     comfortable = None
#     if discomfort_index < -40:
#         discomfort_index_category = 0
#         comfortable = False
#     elif -40 <= discomfort_index < -20:
#         discomfort_index_category = 1
#         comfortable = False
#     elif -20 <= discomfort_index < -10:
#         discomfort_index_category = 2
#         comfortable = False
#     elif -10 <= discomfort_index < -1.8:
#         discomfort_index_category = 3
#         comfortable = False
#     elif -1.8 <= discomfort_index < 13:
#         discomfort_index_category = 4
#         comfortable = False
#     elif 13 <= discomfort_index < 15:
#         discomfort_index_category = 5
#         comfortable = False
#     elif 15 <= discomfort_index < 20:
#         discomfort_index_category = 6
#         comfortable = True
#     elif 20 <= discomfort_index < 26.5:
#         discomfort_index_category = 7
#         comfortable = False
#     elif 26.5 <= discomfort_index < 30:
#         discomfort_index_category = 8
#         comfortable = False
#     elif discomfort_index >= 30:
#         discomfort_index_category = 9
#         comfortable = False
# 
#     return discomfort_index, discomfort_index_category, comfortable
# 
# 
# def wind_chill_temperature(Ta, ws):
#     """
#     Compute the Wind Chill Index, based on the Joint Action Group for Temperature Indices (JAG/TI) methodology
#     :param Ta: temperature in Celcius
#     :type Ta
#     :param ws: wind speed (m/s)
#     :type ws
#     :return: wind chill index, wind chill index category (0-6), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     ws_km_h = ws * 3.6  # convert m/s to km/h wind speed
# 
#     wind_chill_temperature = 13.12 + 0.6215 * Ta - 11.37 * (ws_km_h ** 0.16) + 0.3965 * Ta * (ws_km_h ** 0.16)
# 
#     wind_chill_index_categories = {
#         0: "No risk",
#         1: "Low risk: Slight increase in discomfort",
#         2: "Moderate risk: Risk of hypothermia and frostbite",
#         3: "High risk: High risk of frostbite and hypothermia. Exposed skin can freeze in 10 to 30 minutes",
#         4: "Very high risk: Very high risk of frostbite and hypothermia. Exposed skin can freeze in 5 to 10 minutes",
#         5: "Severe risk: Severe risk of frostbite and hypothermia. Exposed skin can freeze in 2 to 5 minutes",
#         6: "Extreme risk: Outdoor conditions are hazardous. Exposed skin can freeze in less than 2 minutes"
#     }
# 
#     wind_chill_index_category = None
#     comfortable = None
# 
#     if wind_chill_temperature >= 0:
#         wind_chill_index_category = 0
#         comfortable = True
#     elif 0 > wind_chill_temperature >= -9:
#         wind_chill_index_category = 1
#         comfortable = False
#     elif -9 > wind_chill_temperature >= -27:
#         wind_chill_index_category = 2
#         comfortable = False
#     elif -27 > wind_chill_temperature >= -39:
#         wind_chill_index_category = 3
#         comfortable = False
#     elif -39 > wind_chill_temperature >= -47:
#         wind_chill_index_category = 4
#         comfortable = False
#     elif -47 > wind_chill_temperature >= -54:
#         wind_chill_index_category = 5
#         comfortable = False
#     elif wind_chill_temperature < -54:
#         wind_chill_index_category = 6
#         comfortable = False
# 
#     return wind_chill_temperature, wind_chill_index_category, comfortable
# 
# 
# def wet_bulb_globe_temperature_outdoors(Ta, Tr, vp, ws):
#     """
#     Compute the wet-bulb globe temperature from air temperature, mean radiant temperature and vapour pressure
#     :param Ta: temperature in Celcius
#     :type Ta
#     :param Tr: mean radiant temperature (Celsius)
#     :type Tr
#     :param vp: vapour pressure (Pa)
#     :type vp
#     :param ws: wind speed (m/s)
#     :type ws
#     :return: WBGT, WBGT category (0-4), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     wet_bulb_temperature = 1.885 + 0.3704 * Ta + 0.4492 * (vp / 1000)
#     globe_temperature = 2.098 - 2.561 * ws + 0.5957 * Ta + 0.4017 * Tr
#     wet_bulb_globe_temperature = 0.7 * wet_bulb_temperature + 0.2 * globe_temperature + 0.1 * Ta
# 
#     wet_bulb_globe_temperature_categories = {
#         0: "No heat stress",
#         1: "Minor heat stress",
#         2: "Heat stress",
#         3: "High risk: High risk of frostbite and hypothermia. Exposed skin can freeze in 10 to 30 minutes",
#         4: "Very high risk: Very high risk of frostbite and hypothermia. Exposed skin can freeze in 5 to 10 minutes",
#         5: "Severe risk: Severe risk of frostbite and hypothermia. Exposed skin can freeze in 2 to 5 minutes"
#     }
# 
#     wet_bulb_globe_temperature_category = None
#     comfortable = None
# 
#     if wet_bulb_globe_temperature < 18:
#         wet_bulb_globe_temperature_category = 0
#         comfortable = True
#     elif 18 <= wet_bulb_globe_temperature < 23:
#         wet_bulb_globe_temperature_category = 1
#         comfortable = False
#     elif 23 <= wet_bulb_globe_temperature < 28:
#         wet_bulb_globe_temperature_category = 2
#         comfortable = False
#     elif 28 <= wet_bulb_globe_temperature < 30:
#         wet_bulb_globe_temperature_category = 3
#         comfortable = False
#     elif wet_bulb_globe_temperature >= 30:
#         wet_bulb_globe_temperature_category = 4
#         comfortable = False
# 
#     return wet_bulb_globe_temperature, wet_bulb_globe_temperature_category, comfortable
# 
# 
# def effective_temperature(Ta, rh, ws, SR, ac):
#     """
#     Compute the effective temperature
#     :param Ta: temperature (C)
#     :type Ta
#     :param ac: clothing insulation (%)
#     :type ac
#     :param SR: solar radiation (Wh/m2)
#     :type SR
#     :param rh: relative humidity (%)
#     :type rh
#     :param ws: wind speed (m/s)
#     :type ws
#     :return: effective temperature, effective temperature category (0-6), comfortable
#     :rtype: Tuple[float, int, bool]
#     """
# 
#     effective_temperature = None
# 
#     if ws <= 0.2:
#         effective_temperature = Ta - 0.4 * (Ta - 10) * (1 - rh / 100)  # formula by Missenard
#     elif ws > 0.2:
#         effective_temperature = 37 - ((37 - Ta) / (0.68 - (0.0014 * rh) + (1 / (1.76 + 1.4 * (ws ** 0.75))))) - (
#                 0.29 * Ta * (1 - 0.01 * rh))  # modified formula by Gregorczuk (WMO, 1972; Hentschel, 1987)
# 
#     # Radiative-effective temperature
#     radiative_effective_temperature = effective_temperature + ((1 - 0.01 * ac) * SR) * (
#             (0.0155 - 0.00025 * effective_temperature) - (0.0043 - 0.00011 * effective_temperature))
# 
#     effective_temperature_category = None
#     comfortable = None
# 
#     if radiative_effective_temperature < 1:
#         effective_temperature_category = 0
#         comfortable = False
#     elif radiative_effective_temperature >= 1 and effective_temperature < 9:
#         effective_temperature_category = 1
#         comfortable = False
#     elif radiative_effective_temperature >= 9 and effective_temperature < 17:
#         effective_temperature_category = 2
#         comfortable = False
#     elif radiative_effective_temperature >= 17 and effective_temperature < 21:
#         effective_temperature_category = 3
#         comfortable = False
#     elif radiative_effective_temperature >= 21 and effective_temperature < 23:
#         effective_temperature_category = 4
#         comfortable = True
#     elif radiative_effective_temperature >= 23 and effective_temperature < 27:
#         effective_temperature_category = 5
#         comfortable = False
#     elif radiative_effective_temperature >= 27:
#         effective_temperature_category = 6
#         comfortable = False
# 
#     return radiative_effective_temperature, effective_temperature_category, comfortable
# 
# 
# def operative_temperature(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity):
#     import math
#     return (mean_radiant_temperature + (air_temperature * math.sqrt(10 * air_velocity))) / ( 1 + math.sqrt(10 * air_velocity))
# 
# 
# def universal_thermal_climate_index(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity):
#     """
#     Return the approximate Universal Thermal Climate Index for the input criteria.
#     :param air_temperature: Air temperature (C)
#     :param mean_radiant_temperature: Mean radiant temperature (C)
#     :param air_velocity: Air velocity (m/s)
#     :param relative_humidity: Relative humidity
#     :return utci_approx: Approximate UTCI temperature (C)
#     """
#     import math
# 
#     def es(ta):
#         tk = ta + 273.15  # air temp in K
#         es = 2.7150305 * math.log(tk)
#         for count, i in enumerate(
#                 [-2836.5744, -6028.076559, 19.54263612, -0.02737830188, 0.000016261698, (7.0229056 * (10 ** (-10))),
#                  (-1.8680009 * (10 ** (-13)))]):
#             es = es + (i * (tk ** (count - 2)))
#         es = math.exp(es) * 0.01  # convert Pa to hPa
#         return es
# 
#     air_velocity = 0.5 if air_velocity < 0.5 else 17 if air_velocity > 17 else air_velocity
# 
#     # This is a python version of the UTCI_approximation function, Version a 0.002, October 2009
# 
#     ehPa = es(air_temperature) * (relative_humidity / 100.0)
#     D_Tmrt = mean_radiant_temperature - air_temperature
#     Pa = ehPa / 10.0  # convert vapour pressure to kPa
# 
#     utci_approx = air_temperature + \
#                   (0.607562052) + \
#                   (-0.0227712343) * air_temperature + \
#                   (8.06470249 * (10 ** (-4))) * (air_temperature ** 2) + \
#                   (-1.54271372 * (10 ** (-4))) * (air_temperature ** 3) + \
#                   (-3.24651735 * (10 ** (-6))) * (air_temperature ** 4) + \
#                   (7.32602852 * (10 ** (-8))) * (air_temperature ** 5) + \
#                   (1.35959073 * (10 ** (-9))) * (air_temperature ** 6) + \
#                   (-2.25836520) * air_velocity + \
#                   (0.0880326035) * air_temperature * air_velocity + \
#                   (0.00216844454) * (air_temperature ** 2) * air_velocity + \
#                   (-1.53347087 * (10 ** (-5))) * (air_temperature ** 3) * air_velocity + \
#                   (-5.72983704 * (10 ** (-7))) * (air_temperature ** 4) * air_velocity + \
#                   (-2.55090145 * (10 ** (-9))) * (air_temperature ** 5) * air_velocity + \
#                   (-0.751269505) * (air_velocity ** 2) + \
#                   (-0.00408350271) * air_temperature * (air_velocity ** 2) + \
#                   (-5.21670675 * (10 ** (-5))) * (air_temperature ** 2) * (air_velocity ** 2) + \
#                   (1.94544667 * (10 ** (-6))) * (air_temperature ** 3) * (air_velocity ** 2) + \
#                   (1.14099531 * (10 ** (-8))) * (air_temperature ** 4) * (air_velocity ** 2) + \
#                   (0.158137256) * (air_velocity ** 3) + \
#                   (-6.57263143 * (10 ** (-5))) * air_temperature * (air_velocity ** 3) + \
#                   (2.22697524 * (10 ** (-7))) * (air_temperature ** 2) * (air_velocity ** 3) + \
#                   (-4.16117031 * (10 ** (-8))) * (air_temperature ** 3) * (air_velocity ** 3) + \
#                   (-0.0127762753) * (air_velocity ** 4) + \
#                   (9.66891875 * (10 ** (-6))) * air_temperature * (air_velocity ** 4) + \
#                   (2.52785852 * (10 ** (-9))) * (air_temperature ** 2) * (air_velocity ** 4) + \
#                   (4.56306672 * (10 ** (-4))) * (air_velocity ** 5) + \
#                   (-1.74202546 * (10 ** (-7))) * air_temperature * (air_velocity ** 5) + \
#                   (-5.91491269 * (10 ** (-6))) * (air_velocity ** 6) + \
#                   (0.398374029) * D_Tmrt + \
#                   (1.83945314 * (10 ** (-4))) * air_temperature * D_Tmrt + \
#                   (-1.73754510 * (10 ** (-4))) * (air_temperature ** 2) * D_Tmrt + \
#                   (-7.60781159 * (10 ** (-7))) * (air_temperature ** 3) * D_Tmrt + \
#                   (3.77830287 * (10 ** (-8))) * (air_temperature ** 4) * D_Tmrt + \
#                   (5.43079673 * (10 ** (-10))) * (air_temperature ** 5) * D_Tmrt + \
#                   (-0.0200518269) * air_velocity * D_Tmrt + \
#                   (8.92859837 * (10 ** (-4))) * air_temperature * air_velocity * D_Tmrt + \
#                   (3.45433048 * (10 ** (-6))) * (air_temperature ** 2) * air_velocity * D_Tmrt + \
#                   (-3.77925774 * (10 ** (-7))) * (air_temperature ** 3) * air_velocity * D_Tmrt + \
#                   (-1.69699377 * (10 ** (-9))) * (air_temperature ** 4) * air_velocity * D_Tmrt + \
#                   (1.69992415 * (10 ** (-4))) * (air_velocity ** 2) * D_Tmrt + \
#                   (-4.99204314 * (10 ** (-5))) * air_temperature * (air_velocity ** 2) * D_Tmrt + \
#                   (2.47417178 * (10 ** (-7))) * (air_temperature ** 2) * (air_velocity ** 2) * D_Tmrt + \
#                   (1.07596466 * (10 ** (-8))) * (air_temperature ** 3) * (air_velocity ** 2) * D_Tmrt + \
#                   (8.49242932 * (10 ** (-5))) * (air_velocity ** 3) * D_Tmrt + \
#                   (1.35191328 * (10 ** (-6))) * air_temperature * (air_velocity ** 3) * D_Tmrt + \
#                   (-6.21531254 * (10 ** (-9))) * (air_temperature ** 2) * (air_velocity ** 3) * D_Tmrt + \
#                   (-4.99410301 * (10 ** (-6))) * (air_velocity ** 4) * D_Tmrt + \
#                   (-1.89489258 * (10 ** (-8))) * air_temperature * (air_velocity ** 4) * D_Tmrt + \
#                   (8.15300114 * (10 ** (-8))) * (air_velocity ** 5) * D_Tmrt + \
#                   (7.55043090 * (10 ** (-4))) * (D_Tmrt ** 2) + \
#                   (-5.65095215 * (10 ** (-5))) * air_temperature * (D_Tmrt ** 2) + \
#                   (-4.52166564 * (10 ** (-7))) * (air_temperature ** 2) * (D_Tmrt ** 2) + \
#                   (2.46688878 * (10 ** (-8))) * (air_temperature ** 3) * (D_Tmrt ** 2) + \
#                   (2.42674348 * (10 ** (-10))) * (air_temperature ** 4) * (D_Tmrt ** 2) + \
#                   (1.54547250 * (10 ** (-4))) * air_velocity * (D_Tmrt ** 2) + \
#                   (5.24110970 * (10 ** (-6))) * air_temperature * air_velocity * (D_Tmrt ** 2) + \
#                   (-8.75874982 * (10 ** (-8))) * (air_temperature ** 2) * air_velocity * (D_Tmrt ** 2) + \
#                   (-1.50743064 * (10 ** (-9))) * (air_temperature ** 3) * air_velocity * (D_Tmrt ** 2) + \
#                   (-1.56236307 * (10 ** (-5))) * (air_velocity ** 2) * (D_Tmrt ** 2) + \
#                   (-1.33895614 * (10 ** (-7))) * air_temperature * (air_velocity ** 2) * (D_Tmrt ** 2) + \
#                   (2.49709824 * (10 ** (-9))) * (air_temperature ** 2) * (air_velocity ** 2) * (D_Tmrt ** 2) + \
#                   (6.51711721 * (10 ** (-7))) * (air_velocity ** 3) * (D_Tmrt ** 2) + \
#                   (1.94960053 * (10 ** (-9))) * air_temperature * (air_velocity ** 3) * (D_Tmrt ** 2) + \
#                   (-1.00361113 * (10 ** (-8))) * (air_velocity ** 4) * (D_Tmrt ** 2) + \
#                   (-1.21206673 * (10 ** (-5))) * (D_Tmrt ** 3) + \
#                   (-2.18203660 * (10 ** (-7))) * air_temperature * (D_Tmrt ** 3) + \
#                   (7.51269482 * (10 ** (-9))) * (air_temperature ** 2) * (D_Tmrt ** 3) + \
#                   (9.79063848 * (10 ** (-11))) * (air_temperature ** 3) * (D_Tmrt ** 3) + \
#                   (1.25006734 * (10 ** (-6))) * air_velocity * (D_Tmrt ** 3) + \
#                   (-1.81584736 * (10 ** (-9))) * air_temperature * air_velocity * (D_Tmrt ** 3) + \
#                   (-3.52197671 * (10 ** (-10))) * (air_temperature ** 2) * air_velocity * (D_Tmrt ** 3) + \
#                   (-3.36514630 * (10 ** (-8))) * (air_velocity ** 2) * (D_Tmrt ** 3) + \
#                   (1.35908359 * (10 ** (-10))) * air_temperature * (air_velocity ** 2) * (D_Tmrt ** 3) + \
#                   (4.17032620 * (10 ** (-10))) * (air_velocity ** 3) * (D_Tmrt ** 3) + \
#                   (-1.30369025 * (10 ** (-9))) * (D_Tmrt ** 4) + \
#                   (4.13908461 * (10 ** (-10))) * air_temperature * (D_Tmrt ** 4) + \
#                   (9.22652254 * (10 ** (-12))) * (air_temperature ** 2) * (D_Tmrt ** 4) + \
#                   (-5.08220384 * (10 ** (-9))) * air_velocity * (D_Tmrt ** 4) + \
#                   (-2.24730961 * (10 ** (-11))) * air_temperature * air_velocity * (D_Tmrt ** 4) + \
#                   (1.17139133 * (10 ** (-10))) * (air_velocity ** 2) * (D_Tmrt ** 4) + \
#                   (6.62154879 * (10 ** (-10))) * (D_Tmrt ** 5) + \
#                   (4.03863260 * (10 ** (-13))) * air_temperature * (D_Tmrt ** 5) + \
#                   (1.95087203 * (10 ** (-12))) * air_velocity * (D_Tmrt ** 5) + \
#                   (-4.73602469 * (10 ** (-12))) * (D_Tmrt ** 6) + \
#                   (5.12733497) * Pa + \
#                   (-0.312788561) * air_temperature * Pa + \
#                   (-0.0196701861) * (air_temperature ** 2) * Pa + \
#                   (9.99690870 * (10 ** (-4))) * (air_temperature ** 3) * Pa + \
#                   (9.51738512 * (10 ** (-6))) * (air_temperature ** 4) * Pa + \
#                   (-4.66426341 * (10 ** (-7))) * (air_temperature ** 5) * Pa + \
#                   (0.548050612) * air_velocity * Pa + \
#                   (-0.00330552823) * air_temperature * air_velocity * Pa + \
#                   (-0.00164119440) * (air_temperature ** 2) * air_velocity * Pa + \
#                   (-5.16670694 * (10 ** (-6))) * (air_temperature ** 3) * air_velocity * Pa + \
#                   (9.52692432 * (10 ** (-7))) * (air_temperature ** 4) * air_velocity * Pa + \
#                   (-0.0429223622) * (air_velocity ** 2) * Pa + \
#                   (0.00500845667) * air_temperature * (air_velocity ** 2) * Pa + \
#                   (1.00601257 * (10 ** (-6))) * (air_temperature ** 2) * (air_velocity ** 2) * Pa + \
#                   (-1.81748644 * (10 ** (-6))) * (air_temperature ** 3) * (air_velocity ** 2) * Pa + \
#                   (-1.25813502 * (10 ** (-3))) * (air_velocity ** 3) * Pa + \
#                   (-1.79330391 * (10 ** (-4))) * air_temperature * (air_velocity ** 3) * Pa + \
#                   (2.34994441 * (10 ** (-6))) * (air_temperature ** 2) * (air_velocity ** 3) * Pa + \
#                   (1.29735808 * (10 ** (-4))) * (air_velocity ** 4) * Pa + \
#                   (1.29064870 * (10 ** (-6))) * air_temperature * (air_velocity ** 4) * Pa + \
#                   (-2.28558686 * (10 ** (-6))) * (air_velocity ** 5) * Pa + \
#                   (-0.0369476348) * D_Tmrt * Pa + \
#                   (0.00162325322) * air_temperature * D_Tmrt * Pa + \
#                   (-3.14279680 * (10 ** (-5))) * (air_temperature ** 2) * D_Tmrt * Pa + \
#                   (2.59835559 * (10 ** (-6))) * (air_temperature ** 3) * D_Tmrt * Pa + \
#                   (-4.77136523 * (10 ** (-8))) * (air_temperature ** 4) * D_Tmrt * Pa + \
#                   (8.64203390 * (10 ** (-3))) * air_velocity * D_Tmrt * Pa + \
#                   (-6.87405181 * (10 ** (-4))) * air_temperature * air_velocity * D_Tmrt * Pa + \
#                   (-9.13863872 * (10 ** (-6))) * (air_temperature ** 2) * air_velocity * D_Tmrt * Pa + \
#                   (5.15916806 * (10 ** (-7))) * (air_temperature ** 3) * air_velocity * D_Tmrt * Pa + \
#                   (-3.59217476 * (10 ** (-5))) * (air_velocity ** 2) * D_Tmrt * Pa + \
#                   (3.28696511 * (10 ** (-5))) * air_temperature * (air_velocity ** 2) * D_Tmrt * Pa + \
#                   (-7.10542454 * (10 ** (-7))) * (air_temperature ** 2) * (air_velocity ** 2) * D_Tmrt * Pa + \
#                   (-1.24382300 * (10 ** (-5))) * (air_velocity ** 3) * D_Tmrt * Pa + \
#                   (-7.38584400 * (10 ** (-9))) * air_temperature * (air_velocity ** 3) * D_Tmrt * Pa + \
#                   (2.20609296 * (10 ** (-7))) * (air_velocity ** 4) * D_Tmrt * Pa + \
#                   (-7.32469180 * (10 ** (-4))) * (D_Tmrt ** 2) * Pa + \
#                   (-1.87381964 * (10 ** (-5))) * air_temperature * (D_Tmrt ** 2) * Pa + \
#                   (4.80925239 * (10 ** (-6))) * (air_temperature ** 2) * (D_Tmrt ** 2) * Pa + \
#                   (-8.75492040 * (10 ** (-8))) * (air_temperature ** 3) * (D_Tmrt ** 2) * Pa + \
#                   (2.77862930 * (10 ** (-5))) * air_velocity * (D_Tmrt ** 2) * Pa + \
#                   (-5.06004592 * (10 ** (-6))) * air_temperature * air_velocity * (D_Tmrt ** 2) * Pa + \
#                   (1.14325367 * (10 ** (-7))) * (air_temperature ** 2) * air_velocity * (D_Tmrt ** 2) * Pa + \
#                   (2.53016723 * (10 ** (-6))) * (air_velocity ** 2) * (D_Tmrt ** 2) * Pa + \
#                   (-1.72857035 * (10 ** (-8))) * air_temperature * (air_velocity ** 2) * (D_Tmrt ** 2) * Pa + \
#                   (-3.95079398 * (10 ** (-8))) * (air_velocity ** 3) * (D_Tmrt ** 2) * Pa + \
#                   (-3.59413173 * (10 ** (-7))) * (D_Tmrt ** 3) * Pa + \
#                   (7.04388046 * (10 ** (-7))) * air_temperature * (D_Tmrt ** 3) * Pa + \
#                   (-1.89309167 * (10 ** (-8))) * (air_temperature ** 2) * (D_Tmrt ** 3) * Pa + \
#                   (-4.79768731 * (10 ** (-7))) * air_velocity * (D_Tmrt ** 3) * Pa + \
#                   (7.96079978 * (10 ** (-9))) * air_temperature * air_velocity * (D_Tmrt ** 3) * Pa + \
#                   (1.62897058 * (10 ** (-9))) * (air_velocity ** 2) * (D_Tmrt ** 3) * Pa + \
#                   (3.94367674 * (10 ** (-8))) * (D_Tmrt ** 4) * Pa + \
#                   (-1.18566247 * (10 ** (-9))) * air_temperature * (D_Tmrt ** 4) * Pa + \
#                   (3.34678041 * (10 ** (-10))) * air_velocity * (D_Tmrt ** 4) * Pa + \
#                   (-1.15606447 * (10 ** (-10))) * (D_Tmrt ** 5) * Pa + \
#                   (-2.80626406) * (Pa ** 2) + \
#                   (0.548712484) * air_temperature * (Pa ** 2) + \
#                   (-0.00399428410) * (air_temperature ** 2) * (Pa ** 2) + \
#                   (-9.54009191 * (10 ** (-4))) * (air_temperature ** 3) * (Pa ** 2) + \
#                   (1.93090978 * (10 ** (-5))) * (air_temperature ** 4) * (Pa ** 2) + \
#                   (-0.308806365) * air_velocity * (Pa ** 2) + \
#                   (0.0116952364) * air_temperature * air_velocity * (Pa ** 2) + \
#                   (4.95271903 * (10 ** (-4))) * (air_temperature ** 2) * air_velocity * (Pa ** 2) + \
#                   (-1.90710882 * (10 ** (-5))) * (air_temperature ** 3) * air_velocity * (Pa ** 2) + \
#                   (0.00210787756) * (air_velocity ** 2) * (Pa ** 2) + \
#                   (-6.98445738 * (10 ** (-4))) * air_temperature * (air_velocity ** 2) * (Pa ** 2) + \
#                   (2.30109073 * (10 ** (-5))) * (air_temperature ** 2) * (air_velocity ** 2) * (Pa ** 2) + \
#                   (4.17856590 * (10 ** (-4))) * (air_velocity ** 3) * (Pa ** 2) + \
#                   (-1.27043871 * (10 ** (-5))) * air_temperature * (air_velocity ** 3) * (Pa ** 2) + \
#                   (-3.04620472 * (10 ** (-6))) * (air_velocity ** 4) * (Pa ** 2) + \
#                   (0.0514507424) * D_Tmrt * (Pa ** 2) + \
#                   (-0.00432510997) * air_temperature * D_Tmrt * (Pa ** 2) + \
#                   (8.99281156 * (10 ** (-5))) * (air_temperature ** 2) * D_Tmrt * (Pa ** 2) + \
#                   (-7.14663943 * (10 ** (-7))) * (air_temperature ** 3) * D_Tmrt * (Pa ** 2) + \
#                   (-2.66016305 * (10 ** (-4))) * air_velocity * D_Tmrt * (Pa ** 2) + \
#                   (2.63789586 * (10 ** (-4))) * air_temperature * air_velocity * D_Tmrt * (Pa ** 2) + \
#                   (-7.01199003 * (10 ** (-6))) * (air_temperature ** 2) * air_velocity * D_Tmrt * (Pa ** 2) + \
#                   (-1.06823306 * (10 ** (-4))) * (air_velocity ** 2) * D_Tmrt * (Pa ** 2) + \
#                   (3.61341136 * (10 ** (-6))) * air_temperature * (air_velocity ** 2) * D_Tmrt * (Pa ** 2) + \
#                   (2.29748967 * (10 ** (-7))) * (air_velocity ** 3) * D_Tmrt * (Pa ** 2) + \
#                   (3.04788893 * (10 ** (-4))) * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (-6.42070836 * (10 ** (-5))) * air_temperature * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (1.16257971 * (10 ** (-6))) * (air_temperature ** 2) * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (7.68023384 * (10 ** (-6))) * air_velocity * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (-5.47446896 * (10 ** (-7))) * air_temperature * air_velocity * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (-3.59937910 * (10 ** (-8))) * (air_velocity ** 2) * (D_Tmrt ** 2) * (Pa ** 2) + \
#                   (-4.36497725 * (10 ** (-6))) * (D_Tmrt ** 3) * (Pa ** 2) + \
#                   (1.68737969 * (10 ** (-7))) * air_temperature * (D_Tmrt ** 3) * (Pa ** 2) + \
#                   (2.67489271 * (10 ** (-8))) * air_velocity * (D_Tmrt ** 3) * (Pa ** 2) + \
#                   (3.23926897 * (10 ** (-9))) * (D_Tmrt ** 4) * (Pa ** 2) + \
#                   (-0.0353874123) * (Pa ** 3) + \
#                   (-0.221201190) * air_temperature * (Pa ** 3) + \
#                   (0.0155126038) * (air_temperature ** 2) * (Pa ** 3) + \
#                   (-2.63917279 * (10 ** (-4))) * (air_temperature ** 3) * (Pa ** 3) + \
#                   (0.0453433455) * air_velocity * (Pa ** 3) + \
#                   (-0.00432943862) * air_temperature * air_velocity * (Pa ** 3) + \
#                   (1.45389826 * (10 ** (-4))) * (air_temperature ** 2) * air_velocity * (Pa ** 3) + \
#                   (2.17508610 * (10 ** (-4))) * (air_velocity ** 2) * (Pa ** 3) + \
#                   (-6.66724702 * (10 ** (-5))) * air_temperature * (air_velocity ** 2) * (Pa ** 3) + \
#                   (3.33217140 * (10 ** (-5))) * (air_velocity ** 3) * (Pa ** 3) + \
#                   (-0.00226921615) * D_Tmrt * (Pa ** 3) + \
#                   (3.80261982 * (10 ** (-4))) * air_temperature * D_Tmrt * (Pa ** 3) + \
#                   (-5.45314314 * (10 ** (-9))) * (air_temperature ** 2) * D_Tmrt * (Pa ** 3) + \
#                   (-7.96355448 * (10 ** (-4))) * air_velocity * D_Tmrt * (Pa ** 3) + \
#                   (2.53458034 * (10 ** (-5))) * air_temperature * air_velocity * D_Tmrt * (Pa ** 3) + \
#                   (-6.31223658 * (10 ** (-6))) * (air_velocity ** 2) * D_Tmrt * (Pa ** 3) + \
#                   (3.02122035 * (10 ** (-4))) * (D_Tmrt ** 2) * (Pa ** 3) + \
#                   (-4.77403547 * (10 ** (-6))) * air_temperature * (D_Tmrt ** 2) * (Pa ** 3) + \
#                   (1.73825715 * (10 ** (-6))) * air_velocity * (D_Tmrt ** 2) * (Pa ** 3) + \
#                   (-4.09087898 * (10 ** (-7))) * (D_Tmrt ** 3) * (Pa ** 3) + \
#                   (0.614155345) * (Pa ** 4) + \
#                   (-0.0616755931) * air_temperature * (Pa ** 4) + \
#                   (0.00133374846) * (air_temperature ** 2) * (Pa ** 4) + \
#                   (0.00355375387) * air_velocity * (Pa ** 4) + \
#                   (-5.13027851 * (10 ** (-4))) * air_temperature * air_velocity * (Pa ** 4) + \
#                   (1.02449757 * (10 ** (-4))) * (air_velocity ** 2) * (Pa ** 4) + \
#                   (-0.00148526421) * D_Tmrt * (Pa ** 4) + \
#                   (-4.11469183 * (10 ** (-5))) * air_temperature * D_Tmrt * (Pa ** 4) + \
#                   (-6.80434415 * (10 ** (-6))) * air_velocity * D_Tmrt * (Pa ** 4) + \
#                   (-9.77675906 * (10 ** (-6))) * (D_Tmrt ** 2) * (Pa ** 4) + \
#                   (0.0882773108) * (Pa ** 5) + \
#                   (-0.00301859306) * air_temperature * (Pa ** 5) + \
#                   (0.00104452989) * air_velocity * (Pa ** 5) + \
#                   (2.47090539 * (10 ** (-4))) * D_Tmrt * (Pa ** 5) + \
#                   (0.00148348065) * (Pa ** 6)
# 
#     return utci_approx
# 
# 
# def standard_effective_temperature(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity,
#                                    metabolic_rate=1, clo_value=1):
#     """
#     Compute the standard effective temperature
#     :param air_temperature: temperature (C)
#     :type air_temperature
#     :param mean_radiant_temperature: mean radiant temperature (C)
#     :type mean_radiant_temperature
#     :param air_velocity: wind speed (m/s)
#     :type air_velocity
#     :param relative_humidity: relative humidity (%)
#     :type relative_humidity
#     :param metabolic_rate: metabolic rate (met)
#     :type metabolic_rate
#     :param clo_value: clo value ()
#     :type clo_value
#     :return: standard effective temperature
#     :rtype
#     """
#     # TODO: Add source documentation
#     # TODO: Reformat docstring to NumPy style (https://docs.scipy.org/doc/numpy-1.15.0/docs/howto_document.html)
#     # TODO: Add metabolic rates to docstring
# 
#     import math
# 
#     metabolic_rates = {
#         "Sleeping": 0.7,
#         "Reclining": 0.8,
#         "Sitting": 1.0,
#         "Typing": 1.1,
#         "Standing": 1.2,
#         "Driving": 1.5,
#         "Cooking": 1.8,
#         "Walking": 1.7,
#         "Walking 2mph": 2.0,
#         "Lifting 10lbs": 2.1,
#         "Walking 3mph": 2.6,
#         "House Cleaning": 2.7,
#         "Basketball": 3,
#         "Dancing": 3.4,
#         "Walking 4mph": 3.8,
#         "Lifting 100lbs": 4.0,
#         "Shoveling": 4.4,
#         "Running 9mph": 9.5
#     }
# 
#     wme = 0  # External work done - usually around 0 - used for adding additional heat to person given uncontrollable factors
# 
#     # Function to find the saturation vapour pressure, used frequently throughout the comfPierceSET function.
#     def findSaturatedVaporPressureTorr(air_temperature):
#         # calculates Saturated Vapor Pressure (Torr) at Temperature air_temperature  (C)
#         return math.exp(18.6686 - 4030.183 / (air_temperature + 235.0))
# 
#     # Key initial variables.
#     VaporPressure = (relative_humidity * findSaturatedVaporPressureTorr(air_temperature)) / 100
#     AirVelocity = max(air_velocity, 0.1)
#     KCLO = 0.25
#     BODYWEIGHT = 69.9
#     BODYSURFACEAREA = 1.8258
#     METFACTOR = 58.2
#     SBC = 0.000000056697  # Stefan-Boltzmann constant (W/m2K4)
#     CSW = 170
#     CDIL = 120
#     CSTR = 0.5
# 
#     TempSkinNeutral = 33.7  # setpoint (neutral) value for Tsk
#     TempCoreNeutral = 36.8  # setpoint value for Tcr
#     TempBodyNeutral = 36.49  # setpoint for Tb (.1*TempSkinNeutral + .9*TempCoreNeutral)
#     SkinBloodFlowNeutral = 6.3  # neutral value for SkinBloodFlow
# 
#     # INITIAL VALUES - start of 1st experiment
#     TempSkin = TempSkinNeutral
#     TempCore = TempCoreNeutral
#     SkinBloodFlow = SkinBloodFlowNeutral
#     MSHIV = 0.0
#     ALFA = 0.1
#     ESK = 0.1 * metabolic_rate
# 
#     # Start new experiment here (for graded experiments)
#     # UNIT CONVERSIONS (from input variables)
# 
#     p = 101325.0 / 1000  # This variable is the pressure of the atmosphere in kPa and was taken from the psychrometrics.js file of the CBE comfort tool.
# 
#     PressureInAtmospheres = p * 0.009869
#     LTIME = 60
#     TIMEH = LTIME / 60.0
#     RCL = 0.155 * clo_value
#     # AdjustICL(RCL, Conditions);  TH: I don't think this is used in the software
# 
#     FACL = 1.0 + 0.15 * clo_value  # % INCREASE IN BODY SURFACE AREA DUE TO CLOTHING
#     LR = 2.2 / PressureInAtmospheres  # Lewis Relation is 2.2 at sea level
#     RM = metabolic_rate * METFACTOR
#     M = metabolic_rate * METFACTOR
# 
#     if clo_value <= 0:
#         WCRIT = 0.38 * pow(AirVelocity, -0.29)
#         ICL = 1.0
#     else:
#         WCRIT = 0.59 * pow(AirVelocity, -0.08)
#         ICL = 0.45
# 
#     CHC = 3.0 * pow(PressureInAtmospheres, 0.53)
#     CHCV = 8.600001 * pow((AirVelocity * PressureInAtmospheres), 0.53)
#     CHC = max(CHC, CHCV)
# 
#     # initial estimate of Tcl
#     CHR = 4.7
#     CTC = CHR + CHC
#     RA = 1.0 / (FACL * CTC)  # resistance of air layer to dry heat transfer
#     TOP = (CHR * mean_radiant_temperature + CHC * air_temperature) / CTC
#     TCL = TOP + (TempSkin - TOP) / (CTC * (RA + RCL))
# 
#     # ========================  BEGIN ITERATION
#     #
#     # Tcl and CHR are solved iteratively using: H(Tsk - To) = CTC(Tcl - To),
#     # where H = 1/(Ra + Rcl) and Ra = 1/Facl*CTC
# 
#     TCL_OLD = TCL
#     TIME = range(LTIME)
#     flag = True
#     for TIM in TIME:
#         if flag == True:
#             while abs(TCL - TCL_OLD) > 0.01:
#                 TCL_OLD = TCL
#                 CHR = 4.0 * SBC * pow(((TCL + mean_radiant_temperature) / 2.0 + 273.15), 3.0) * 0.72
#                 CTC = CHR + CHC
#                 RA = 1.0 / (FACL * CTC)  # resistance of air layer to dry heat transfer
#                 TOP = (CHR * mean_radiant_temperature + CHC * air_temperature) / CTC
#                 TCL = (RA * TempSkin + RCL * TOP) / (RA + RCL)
#         flag = False
#         DRY = (TempSkin - TOP) / (RA + RCL)
#         HFCS = (TempCore - TempSkin) * (5.28 + 1.163 * SkinBloodFlow)
#         ERES = 0.0023 * M * (44.0 - VaporPressure)
#         CRES = 0.0014 * M * (34.0 - air_temperature)
#         SCR = M - HFCS - ERES - CRES - wme
#         SSK = HFCS - DRY - ESK
#         TCSK = 0.97 * ALFA * BODYWEIGHT
#         TCCR = 0.97 * (1 - ALFA) * BODYWEIGHT
#         DTSK = (SSK * BODYSURFACEAREA) / (TCSK * 60.0)  # //deg C per minute
#         DTCR = SCR * BODYSURFACEAREA / (TCCR * 60.0)  # //deg C per minute
#         TempSkin = TempSkin + DTSK
#         TempCore = TempCore + DTCR
#         TB = ALFA * TempSkin + (1 - ALFA) * TempCore
#         SKSIG = TempSkin - TempSkinNeutral
#         WARMS = (SKSIG > 0) * SKSIG
#         COLDS = ((-1.0 * SKSIG) > 0) * (-1.0 * SKSIG)
#         CRSIG = (TempCore - TempCoreNeutral)
#         WARMC = (CRSIG > 0) * CRSIG
#         COLDC = ((-1.0 * CRSIG) > 0) * (-1.0 * CRSIG)
#         BDSIG = TB - TempBodyNeutral
#         WARMB = (BDSIG > 0) * BDSIG
#         COLDB = ((-1.0 * BDSIG) > 0) * (-1.0 * BDSIG)
#         SkinBloodFlow = (SkinBloodFlowNeutral + CDIL * WARMC) / (1 + CSTR * COLDS)
#         if SkinBloodFlow > 90.0: SkinBloodFlow = 90.0
#         if SkinBloodFlow < 0.5: SkinBloodFlow = 0.5
#         REGSW = CSW * WARMB * math.exp(WARMS / 10.7)
#         if REGSW > 500.0: REGSW = 500.0
#         ERSW = 0.68 * REGSW
#         REA = 1.0 / (LR * FACL * CHC)  # evaporative resistance of air layer
#         RECL = RCL / (LR * ICL)  # evaporative resistance of clothing (icl=.45)
#         EMAX = (findSaturatedVaporPressureTorr(TempSkin) - VaporPressure) / (REA + RECL)
#         PRSW = ERSW / EMAX
#         PWET = 0.06 + 0.94 * PRSW
#         EDIF = PWET * EMAX - ERSW
#         ESK = ERSW + EDIF
#         if PWET > WCRIT:
#             PWET = WCRIT
#             PRSW = WCRIT / 0.94
#             ERSW = PRSW * EMAX
#             EDIF = 0.06 * (1.0 - PRSW) * EMAX
#             ESK = ERSW + EDIF
#         if EMAX < 0:
#             EDIF = 0
#             ERSW = 0
#             PWET = WCRIT
#             PRSW = WCRIT
#             ESK = EMAX
#         ESK = ERSW + EDIF
#         MSHIV = 19.4 * COLDS * COLDC
#         M = RM + MSHIV
#         ALFA = 0.0417737 + 0.7451833 / (SkinBloodFlow + .585417)
# 
#     # Define new heat flow terms, coeffs, and abbreviations
#     STORE = M - wme - CRES - ERES - DRY - ESK  # rate of body heat storage
#     HSK = DRY + ESK  # total heat loss from skin
#     RN = M - wme  # net metabolic heat production
#     ECOMF = 0.42 * (RN - (1 * METFACTOR))
#     if ECOMF < 0.0: ECOMF = 0.0  # from Fanger
#     EREQ = RN - ERES - CRES - DRY
#     EMAX = EMAX * WCRIT
#     HD = 1.0 / (RA + RCL)
#     HE = 1.0 / (REA + RECL)
#     W = PWET
#     PSSK = findSaturatedVaporPressureTorr(TempSkin)
#     # Definition of ASHRAE standard environment... denoted "S"
#     CHRS = CHR
#     if metabolic_rate < 0.85:
#         CHCS = 3.0
#     else:
#         CHCS = 5.66 * pow((metabolic_rate - 0.85), 0.39)
#         if CHCS < 3.0: CHCS = 3.0
# 
#     CTCS = CHCS + CHRS
#     RCLOS = 1.52 / ((metabolic_rate - wme / METFACTOR) + 0.6944) - 0.1835
#     RCLS = 0.155 * RCLOS
#     FACLS = 1.0 + KCLO * RCLOS
#     FCLS = 1.0 / (1.0 + 0.155 * FACLS * CTCS * RCLOS)
#     IMS = 0.45
#     ICLS = IMS * CHCS / CTCS * (1 - FCLS) / (CHCS / CTCS - FCLS * IMS)
#     RAS = 1.0 / (FACLS * CTCS)
#     REAS = 1.0 / (LR * FACLS * CHCS)
#     RECLS = RCLS / (LR * ICLS)
#     HD_S = 1.0 / (RAS + RCLS)
#     HE_S = 1.0 / (REAS + RECLS)
# 
#     # SET* (standardized humidity, clo_value, Pb, and CHC) determined using Newton's iterative solution
#     # FNERRS is defined in the GENERAL SETUP section above
# 
#     DELTA = .0001
#     dx = 100.0
#     X_OLD = TempSkin - HSK / HD_S  # lower bound for SET
#     while abs(dx) > .01:
#         ERR1 = (HSK - HD_S * (TempSkin - X_OLD) - W * HE_S * (PSSK - 0.5 * findSaturatedVaporPressureTorr(X_OLD)))
#         ERR2 = (HSK - HD_S * (TempSkin - (X_OLD + DELTA)) - W * HE_S * (
#                 PSSK - 0.5 * findSaturatedVaporPressureTorr((X_OLD + DELTA))))
#         X = X_OLD - DELTA * ERR1 / (ERR2 - ERR1)
#         dx = X - X_OLD
#         X_OLD = X
# 
#     return X
# 
# 
# def pmv_elevated_airspeed(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity, metabolic_rate=1,
#                           clo_value=1):
#     """
#     Compute the Positive Mean Vote (PMV) with an attempt to incorporate the cooling effect of air movement beyond
#     0.15m/s, testing against the Standard Effective Temperature method. If air velocity is below a threshold value, the
#     original Fanger PMV method will be used.
#     :param air_temperature:
#     :param mean_radiant_temperature:
#     :param air_velocity:
#     :param relative_humidity:
#     :param metabolic_rate:
#     :param clo_value:
#     :return pmv: Predicted mean vote
#     :return ppd: Percent predicted dissatisfied [%]
#     """
# 
#     r = list()
#     s_e_t = standard_effective_temperature(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity,
#                                            metabolic_rate, clo_value)
#     still_air_threshold = 0.1
# 
#     def util_secant(a, b, epsilon):
#         """
#         Function taken from the util.js script of the CBE comfort tool page and modified to include the fn inside the
#         util_secant function.
#         :param a:
#         :param b:
#         :param epsilon:
#         :return:
#         """
#         # root-finding only
#         res = []
# 
#         def fn(t):
#             if t > 200:
#                 t = 200
#             if t < -200:
#                 t = 200
#             return s_e_t - standard_effective_temperature(
#                 air_temperature - t, mean_radiant_temperature - t, still_air_threshold, relative_humidity,
#                 metabolic_rate, clo_value)
# 
#         f1 = fn(a)
#         if abs(f1) <= epsilon:
#             res.append(a)
#         else:
#             f2 = fn(b)
#             if abs(f2) <= epsilon:
#                 res.append(b)
#             else:
#                 count = range(100)
#                 for i in count:
#                     if (b - a) != 0 and (f2 - f1) != 0:
#                         slope = (f2 - f1) / (b - a)
#                         c = b - f2 / slope
#                         f3 = fn(c)
#                         if abs(f3) < epsilon:
#                             res.append(c)
#                         a = b
#                         b = c
#                         f1 = f2
#                         f2 = f3
#             res.append('NaN')
# 
#         return res[0]
# 
#     def util_bisect(a, b, epsilon, target):
#         """
#         Function taken from the util.js script of the CBE comfort tool page and modified to include the fn inside the
#         util_secant function definition.
#         :param a:
#         :param b:
#         :param epsilon:
#         :param target:
#         :return:
#         """
# 
#         def fn(t):
#             return s_e_t - standard_effective_temperature(
#                 air_temperature - t, mean_radiant_temperature - t, still_air_threshold, relative_humidity,
#                 metabolic_rate, clo_value)
# 
#         while abs(b - a) > (2 * epsilon):
#             midpoint = (b + a) / 2
#             a_t = fn(a)
#             b_t = fn(b)
#             midpoint_t = fn(midpoint)
#             if (a_t - target) * (midpoint_t - target) < 0:
#                 b = midpoint
#             elif (b_t - target) * (midpoint_t - target) < 0:
#                 a = midpoint
#             else:
#                 return -999
#         return midpoint
# 
#     if air_velocity <= still_air_threshold:
#         pmv, ppd = fanger_pmv(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity,
#                               metabolic_rate, clo_value)
#         ta_adj = air_temperature
#         ce = 0
#     else:
#         ce_l = 0
#         ce_r = 40
#         eps = 0.001  # precision of ce
# 
#         ce = util_secant(ce_l, ce_r, eps)
#         if ce == 'NaN':
#             ce = util_bisect(ce_l, ce_r, eps, 0)
# 
#         pmv, ppd = fanger_pmv(air_temperature - ce, mean_radiant_temperature - ce, still_air_threshold,
#                               relative_humidity, metabolic_rate, clo_value)
#         ta_adj = air_temperature - ce
#         tr_adj = mean_radiant_temperature - ce
# 
#     r.append(pmv)
#     r.append(ppd)
#     # r.append(s_e_t)
#     # r.append(ta_adj)
#     # r.append(ce)
# 
#     return r
# 
# 
# def fanger_pmv(air_temperature, mean_radiant_temperature, air_velocity, relative_humidity, metabolic_rate=1,
#                clo_value=1):
#     """
#     Original Fanger function to compute PMV.  Only intended for use with low air
#     speeds (<0.1 m/s).
#     Args:
#         air_temperature: air temperature (C)
#         mean_radiant_temperature: mean radiant temperature (C)
#         air_velocity: relative air velocity (m/s)
#         relative_humidity: relative humidity (%) Used only this way to input humidity level
#         metabolic_rate: metabolic rate (metabolic_rate)
#         clo_value: clothing (clo)
#     Returns:
#         pmv: predicted mean vote
#         ppd: percentage of people dissatisfied.
#     """
# 
#     pa = relative_humidity * 10 * math.exp(16.6536 - 4030.183 / (air_temperature + 235))
# 
#     icl = 0.155 * clo_value  # thermal insulation of the clothing in M2K/W
#     m = metabolic_rate * 58.15  # metabolic rate in W/M2
#     wme = 0  # Additional work done - no idea why this is here!
#     w = wme * 58.15  # external work in W/M2
#     mw = m - w  # internal heat production in the human body
#     if (icl <= 0.078):
#         fcl = 1 + (1.29 * icl)
#     else:
#         fcl = 1.05 + (0.645 * icl)
# 
#     # heat transf. coeff. by forced convection
#     hcf = 12.1 * math.sqrt(air_velocity)
#     taa = air_temperature + 273
#     tra = mean_radiant_temperature + 273
#     tcla = taa + (35.5 - air_temperature) / (3.5 * icl + 0.1)
# 
#     p1 = icl * fcl
#     p2 = p1 * 3.96
#     p3 = p1 * 100
#     p4 = p1 * taa
#     p5 = (308.7 - 0.028 * mw) + (p2 * math.pow(tra / 100, 4))
#     xn = tcla / 100
#     xf = tcla / 50
#     eps = 0.00015
# 
#     n = 0
#     while abs(xn - xf) > eps:
#         xf = (xf + xn) / 2
#         hcn = 2.38 * math.pow(abs(100.0 * xf - taa), 0.25)
#         if (hcf > hcn):
#             hc = hcf
#         else:
#             hc = hcn
#         xn = (p5 + p4 * hc - p2 * math.pow(xf, 4)) / (100 + p3 * hc)
#         n += 1
#         if (n > 150):
#             print('Max iterations exceeded')
#             return 1
# 
#     tcl = 100 * xn - 273
# 
#     # heat loss diff. through skin
#     hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa)
#     # heat loss by sweating
#     if mw > 58.15:
#         hl2 = 0.42 * (mw - 58.15)
#     else:
#         hl2 = 0
#     # latent respiration heat loss
#     hl3 = 1.7 * 0.00001 * m * (5867 - pa)
#     # dry respiration heat loss
#     hl4 = 0.0014 * m * (34 - air_temperature)
#     # heat loss by radiation
#     hl5 = 3.96 * fcl * (math.pow(xn, 4) - math.pow(tra / 100, 4))
#     # heat loss by convection
#     hl6 = fcl * hc * (tcl - air_temperature)
# 
#     ts = 0.303 * math.exp(-0.036 * m) + 0.028
#     pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)
#     ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(pmv, 4.0) - 0.2179 * pow(pmv, 2.0))
# 
#     r = list()
#     r.append(pmv)
#     r.append(ppd)
# 
#     return r
# 
# 
# def degree_days_heating(external_temperatures, heating_temperature=18.3):
#     """
#     Compute the number of heating degree days for each day based on the given heating reference temperature, and using
#     the most accurate mean degree-hours; calculated from the hourly temperature record method:
#     :param external_temperatures: temperature in Celcius
#     :type external_temperatures: list(float) (n=8760)
#     :param heating_temperature: heating reference temperature
#     :type heating_temperature
#     :return: heating degree days per day
#     :rtype: [float] (n=365)
#     """
# 
#     # Get heating degree hours
#     heating_degree_hours = [heating_temperature - i for i in external_temperatures]
#     # Replace negative values with 0
#     heating_degree_hours_positive = [i if i >= 0 else 0 for i in heating_degree_hours]
#     # Convert to degree-days
#     heating_degree_hours_hourly = [i / 24 for i in heating_degree_hours_positive]
#     # Sum daily values to get degree-days per day
#     heating_degree_days_daily = [sum(i) for i in chunks(heating_degree_hours_hourly, 24)]
# 
#     return heating_degree_days_daily
# 
# 
# def degree_days_cooling(external_temperatures, cooling_temperature=23.3):
#     """
#     Compute the number of cooling degree days for each day based on the given cooling reference temperature, and using
#     the most accurate mean degree-hours; calculated from the hourly temperature record method:
#     :param external_temperatures: temperature in Celcius
#     :type external_temperatures: list(float) (n=8760)
#     :param cooling_temperature: cooling reference temperature
#     :type cooling_temperature
#     :return: heating degree days per day
#     :rtype: [float] (n=365)
#     """
# 
#     # Get heating degree hours
#     cooling_degree_hours = [i - cooling_temperature for i in external_temperatures]
#     # Replace negative values with 0
#     cooling_degree_hours_positive = [i if i >= 0 else 0 for i in cooling_degree_hours]
#     # Convert to degree-days
#     cooling_degree_days_hourly = [i / 24 for i in cooling_degree_hours_positive]
#     # Sum daily values to get degree-days per day
#     cooling_degree_days_daily = [sum(i) for i in chunks(cooling_degree_days_hourly, 24)]
# 
#     return cooling_degree_days_daily
# 
# 
# def celsius_to_kelvin(t):
#     """ Converts degrees Celsius into Kelvin
# 
#     Example
#     >>> celsius_to_kelvin(0)
#     273.15
# 
#     :param t: temperature (C)
#     :return: temperature (K)
#     """
#     return t + 273.15
# 
# 
# def pressure_at_altitude(altitude):
#     """ Calculate pressure as a function of altitude
# 
#     Source: EQ(3) 2009 ASHRAE Handbook—Fundamentals (SI)
# 
#     :param altitude: altitude in m
#     :return: pressure in kPa
#     """
#     return 101.325 * (1 - 0.0000225577 * altitude) ** 5.2559
# 
# 
# def temperature_at_altitude(altitude):
#     """ Calculate temperature as a function of altitude
# 
#     Source: EQ(4) 2009 ASHRAE Handbook—Fundamentals (SI)
# 
#     :param altitude: altitude in m
#     :return: temperature (C)
#     """
#     return 15 - 0.0065 * altitude
# 
# 
# def temperature_enthalpy_humidity_ratio(enthalpy, humidity_ratio):
#     """ Calculate temperature as a function of enthalpy and humidity ratio
# 
#     Source: UNKNOWN CURRENTLY
# 
#     :param enthalpy: enthalpy in m
#     :param humidity_ratio: humidity ratio in g/g
#     :return: temperature (C)
#     """
#     return (enthalpy - 2.5 * (humidity_ratio * 1000)) / (1.01 + (0.00189 * humidity_ratio * 1000))
# 
# 
# def partial_vapour_pressure(dry_bulb_temperature):
#     """ Calculate saturation pressure over liquid water for the temperature range -100C to 200C
# 
#     Source: EQ(5) & EQ(6) 2009 ASHRAE Handbook—Fundamentals (SI)
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :return: Saturation pressure (kPa)
#     """
#     c1 = -5674.5359
#     c2 = 6.3925247
#     c3 = -0.009677843
#     c4 = 0.00000062215701
#     c5 = 2.0747825E-09
#     c6 = -9.484024E-13
#     c7 = 4.1635019
#     c8 = -5800.2206
#     c9 = 1.3914993
#     c10 = -0.048640239
#     c11 = 0.000041764768
#     c12 = -0.000000014452093
#     c13 = 6.5459673
# 
#     t_k = celsius_to_kelvin(dry_bulb_temperature)  # convert Celsius into Kelvin
# 
#     if dry_bulb_temperature < 0:
#         return math.exp(c1 / t_k + c2 + c3 * t_k + c4 * t_k ** 2 + c5 * t_k ** 3 + c6 * t_k ** 4 + c7 * math.log(t_k))
#     else:
#         return math.exp(c8 / t_k + c9 + c10 * t_k + c11 * t_k ** 2 + c12 * t_k ** 3 + c13 * math.log(t_k))
# 
# 
# def density_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure):
#     """ Calculate density from dry bulb temperature and wet bulb temp
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: density (kg/m3)
#     """
#     return (1 / specific_volume_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure)) * (
#                 1 + humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure))
# 
# 
# def specific_volume_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature,
#                                          pressure):
#     """ Calculate Specific Volume from dry bulb temperature and wet bulb temp
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: Specific Volume (m3/kg)
#     """
#     return 287.042 * celsius_to_kelvin(dry_bulb_temperature) * (
#                 1 + 1.607858 * humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature,
#                                                                    pressure)) / (pressure * 1000)
# 
# 
# def humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature,
#                                         pressure):
#     """ Calculate Humidity Ratio from dry bulb temperature & wet bulb temp
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: humidity_ratio (kg_water/kg_dry_air)
#     """
#     # calculate saturation humidity ratio of wet bulb temperature
#     ws = 0.621945 * partial_vapour_pressure(wet_bulb_temperature) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature))
#     return ((2501 - 2.326 * wet_bulb_temperature) * ws - 1.006 * (dry_bulb_temperature - wet_bulb_temperature)) / (
#                 2501 + 1.86 * dry_bulb_temperature - 4.186 * wet_bulb_temperature)
# 
# 
# def relative_humidity_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature,
#                                            pressure):
#     """ Calculate Relative Humidity from dry bulb and wet bulb temperature
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: Relative humidity (%)
#     """
#     # calculate partial vapour pressure
#     pvp = humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure) * pressure * 1000 / (
#                 humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure) + 0.621945)
# 
#     return min([1, pvp / partial_vapour_pressure(dry_bulb_temperature)])
# 
# 
# def dew_point_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure):
#     """ Calculate dew temperature from dry bulb temperature and wet bulb temperature
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: dew point temperature (C)
#     """
#     c14 = 6.54
#     c15 = 14.526
#     c16 = 0.7389
#     c17 = 0.09486
#     c18 = 0.4569
# 
#     wg = humidity_ratio_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure)
# 
#     # calculate partial vapour pressure
#     pvp = pressure * 1000 * wg / (0.621945 + wg)
#     a = math.log(pvp / 1000)
# 
#     if dry_bulb_temperature < 0:
#         return 6.09 + 12.608 * a + 0.4959 * a ** 2
#     else:
#         return c14 + c15 * a + c16 * a ** 2 + c17 * a ** 3 + c18 * (pvp / 1000) ** 0.1984
# 
# 
# def enthalpy_wet_bulb_temperature(dry_bulb_temperature, wet_bulb_temperature, pressure):
#     """ Calculate enthalpy from dry bulb temperature and wet bulb temperature
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param wet_bulb_temperature: wet-bulb temperature (C)
#     :param pressure: pressure (kPa)
#     :return: enthalpy (J/kg)
#     """
#     return 1.006 * dry_bulb_temperature + humidity_ratio_wet_bulb_temperature(dry_bulb_temperature,
#                                                                               wet_bulb_temperature, pressure) * (
#                        2501 + 1.86 * dry_bulb_temperature)
# 
# 
# def density_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
#     """ Calculate density from dry bulb temperature and relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param pressure: pressure (kPa)
#     :return: density (kg/m3)
#     """
#     return (1 / specific_volume_relative_humidity(dry_bulb_temperature, relative_humidity, pressure)) * (
#                 1 + humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure))
# 
# 
# def specific_volume_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
#     """ Calculate Specific Volume from dry bulb temperature and relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param pressure: pressure (kPa)
#     :return: Specific Volume (m3/kg)
#     """
#     return 287.042 * celsius_to_kelvin(dry_bulb_temperature) * (
#                 1 + 1.607858 * humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure)) / (
#                        pressure * 1000)
# 
# 
# def wet_bulb_temperature_relative_humidity(dry_bulb_temperature, relative_humidity,
#                                            pressure):
#     """ Calculate wet bulb temperature from dry bulb temperature and relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param pressure: pressure (kPa)
#     :return: Wet-bulb temperature (C)
#     """
# 
#     # Initial Wet Bulb Temp (equal to dry_bulb_temperature)
#     wet_bulb_temperature1 = dry_bulb_temperature
# 
#     # calculate saturation humidity ratio of wet bulb temperature
#     saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
# 
#     # calculate humidity ratio of partial vapour pressure
#     humidity_ratio = humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure)
# 
#     # calculate wet bulb temperature
#     wet_bulb_temperature2 = ((2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 4.186 * humidity_ratio - 2.326 * saturation_vapour_pressure)
# 
#     myerror = wet_bulb_temperature2 - wet_bulb_temperature1
# 
#     if dry_bulb_temperature < 0:
#         while abs(myerror) > 0.001:
#             wet_bulb_temperature1 = 0.1 * myerror + wet_bulb_temperature1
#             saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
#             wet_bulb_temperature2 = ((2830 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2830 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 2.1 * humidity_ratio - 0.24 * saturation_vapour_pressure)
#             myerror = wet_bulb_temperature2 - wet_bulb_temperature1
#     else:
#         while abs(myerror) > 0.001:
#             wet_bulb_temperature1 = 0.01 * myerror + wet_bulb_temperature1
#             saturation_vapour_pressure = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
#             wet_bulb_temperature2 = ((2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * saturation_vapour_pressure + 1.006 * dry_bulb_temperature) / (1.006 + 4.186 * humidity_ratio - 2.326 * saturation_vapour_pressure)
#             myerror = wet_bulb_temperature2 - wet_bulb_temperature1
# 
#     return wet_bulb_temperature2
# 
# 
# def humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
#     """ Calculate Humidity Ratio from dry bulb temperature & relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param pressure: pressure (kPa)
#     :return: humidity_ratio (kg_water/kg_dry_air)
#     """
#     # calculate partial water pressure
#     pvp = relative_humidity * partial_vapour_pressure(dry_bulb_temperature)
# 
#     # calculate humidity ratio
#     return 0.621945 * pvp / (pressure * 1000 - pvp)
# 
# 
# def dew_point_relative_humidity(dry_bulb_temperature, relative_humidity):
#     """ Calculate dew point temperature from dry bulb temperature and relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :return: dew point temperature (C)
#     """
# 
#     c14 = 6.54
#     c15 = 14.526
#     c16 = 0.7389
#     c17 = 0.09486
#     c18 = 0.4569
# 
#     # calculate partial vapour pressure
#     pvp = relative_humidity * partial_vapour_pressure(dry_bulb_temperature)
#     a = math.log(pvp / 1000)
# 
#     if dry_bulb_temperature < 0:
#         return 6.09 + 12.608 * a + 0.4959 * a ** 2
#     else:
#         return c14 + c15 * a + c16 * a ** 2 + c17 * a ** 3 + c18 * (pvp / 1000) ** 0.1984
# 
# 
# def enthalpy_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
#     """ Calculate enthalpy from dry bulb temperature and relative humidity
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param relative_humidity: Relative humidity (%)
#     :param pressure: pressure (kPa)
#     :return: humidity_ratio (kg_water/kg_dry_air)
#     """
#     return 1.005 * dry_bulb_temperature + humidity_ratio_relative_humidity(dry_bulb_temperature, relative_humidity,
#                                                                            pressure) * ( 2501 + 1.86 * dry_bulb_temperature)
# 
# 
# def density_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
#     """ Calculate density from dry bulb temperature and Humidity Ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :param pressure: pressure (kPa)
#     :return: density (kg/m3)
#     """
#     return (1 / specific_volume_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure)) * (1 + humidity_ratio)
# 
# 
# def specific_volume_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
#     """ Calculate Specific Volume from dry bulb temperature and Humidity Ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :param pressure: pressure (kPa)
#     :return: Specific volume (m3/kg)
#     """
#     return 287.042 * celsius_to_kelvin(dry_bulb_temperature) * (1 + 1.607858 * humidity_ratio) / (pressure * 1000)
# 
# 
# def wet_bulb_temperature_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
#     """ Calculate wet bulb temperature from dry bulb temperature and humidity ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :param pressure: pressure (kPa)
#     :return: Wet-bulb temperature (C)
#     """
# 
#     # Initial wet_bulb_temperature_rh (equal to dbt)
#     wet_bulb_temperature1 = dry_bulb_temperature
# 
#     # calculate saturation humidity ratio of dry bulb temperature
#     ws1 = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
# 
#     # calculate wet bulb temperature
#     wet_bulb_temperature2 = ((2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * ws1 + 1.006 * dry_bulb_temperature) / (1.006 + 4.186 * humidity_ratio - 2.326 * ws1)
# 
#     my_error = wet_bulb_temperature2 - wet_bulb_temperature1
# 
#     if dry_bulb_temperature < 0:
#         while abs(my_error) > 0.001:
#             wet_bulb_temperature1 = 0.1 * my_error + wet_bulb_temperature1
#             ws1 = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
#             wet_bulb_temperature2 = ((
#                                                  2830 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2830 * ws1 + 1.006 * dry_bulb_temperature) / (
#                                                 1.006 + 2.1 * humidity_ratio - 0.24 * ws1)
#             my_error = wet_bulb_temperature2 - wet_bulb_temperature1
#     else:
#         while abs(my_error) > 0.001:
#             wet_bulb_temperature1 = 0.01 * my_error + wet_bulb_temperature1
#             ws1 = 0.621945 * partial_vapour_pressure(wet_bulb_temperature1) / (pressure * 1000 - partial_vapour_pressure(wet_bulb_temperature1))
#             wet_bulb_temperature2 = ((
#                                                  2501 + 1.86 * dry_bulb_temperature) * humidity_ratio - 2501 * ws1 + 1.006 * dry_bulb_temperature) / (
#                                                 1.006 + 4.186 * humidity_ratio - 2.326 * ws1)
#             my_error = wet_bulb_temperature2 - wet_bulb_temperature1
#     return wet_bulb_temperature2
# 
# 
# def relative_humidity_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
#     """ Calculate Relative Humidity from dry bulb temperature and humidity ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :param pressure: pressure (kPa)
#     :return: Relative humidity (%)
#     """
# 
#     # calculate partial vapour pressure
#     pvp = humidity_ratio * pressure * 1000 / (humidity_ratio + 0.621945)
# 
#     return min([1, pvp / partial_vapour_pressure(dry_bulb_temperature)])
# 
# 
# def dew_point_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
#     """ Calculate dew point temperature from dry bulb temperature and humidity ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :param pressure: pressure (kPa)
#     :return: dew point temperature (C)
#     """
# 
#     c14 = 6.54
#     c15 = 14.526
#     c16 = 0.7389
#     c17 = 0.09486
#     c18 = 0.4569
# 
#     # calculate partial vapour pressure
#     pvp = humidity_ratio * pressure * 1000 / (humidity_ratio + 0.621945)
#     a = math.log(pvp / 1000)
# 
#     if dry_bulb_temperature < 0:
#         return 6.09 + 12.608 * a + 0.4959 * a ** 2
#     else:
#         return c14 + c15 * a + c16 * a ** 2 + c17 * a ** 3 + c18 * (pvp / 1000) ** 0.1984
# 
# 
# def enthalpy_humidity_ratio(dry_bulb_temperature, humidity_ratio):
#     """ Calculate enthalpy from dry bulb temperature and humidity ratio
# 
#     :param dry_bulb_temperature: dry-bulb temperature (C)
#     :param humidity_ratio: humidity_ratio (kg_water/kg_dry_air)
#     :return: dew point temperature (C)
#     """
#     return 1.006 * dry_bulb_temperature + humidity_ratio * (2501 + 1.86 * dry_bulb_temperature)