# import numpy as np
# from .constants import KELVIN, STEFAN_BOLTZMANN_CONSTANT
#
#
# def calc_sky_emissivity(dew_point_temperature, total_sky_cover):
#     """
#     Calculate the sky emissivity using the Clark-Allen method used by EnergyPlus.
#
#     Clark and Allen report that the standard estimation error for the atmospheric long-wave radiation was 10 𝑊𝑊/𝑚𝑚2.
#
#     Parameters
#     ----------
#     dew_point_temperature : float
#         Dew-point temperature of air in degrees C
#     total_sky_cover : int
#         Total sky cover in tenths. Conversion to decimal 0-1 happens within the method.
#
#     Returns
#     -------
#     sky_emissivity : float
#
#     """
#
#     sky_emissivity = (0.787 + 0.764 * np.log((dew_point_temperature + KELVIN) / KELVIN)) * (
#             1 + 0.0224 * total_sky_cover - 0.0035 * np.power(total_sky_cover, 2) + 0.0028 * np.power(
#         total_sky_cover, 3))
#     return sky_emissivity
#
#
# def simple_combined_exterior_heat_transfer_coefficient(material_roughness, local_wind_speed):
#     """
#     The simple algorithm uses surface roughness and local surface windspeed to calculate the exterior heat transfer coefficient.
#
#     Note that the simple correlation yields a combined convection and radiation heat transfer coefficient.
#
#     Material roughness correlations are taken from Figure 1, Page 22.4, ASHRAE Handbook of Fundamentals (ASHRAE 1989).
#
#     Source: http://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-020.html for further details.
#
#     Parameters
#     ----------
#     material_roughness : string
#         A surface roughness index. Acceptable values are: "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth" and "Very smooth"
#     local_wind_speed : float
#         Local wind speed calculated at the height above ground of the surface centroid
#
#     Returns
#     -------
#     heat_transfer_coefficient
#         The surface heat transfer coefficient for the variables passed
#     """
#
#     material = {
#         "Very rough": {"D": 11.58, "E": 5.89, "F": 0},
#         "Rough": {"D": 12.49, "E": 4.065, "F": 0.028},
#         'Medium rough': {"D": 10.79, "E": 4.192, "F": 0},
#         'Medium smooth': {"D": 8.23, "E": 4, "F": -0.057},
#         'Smooth': {"D": 10.22, "E": 3.1, "F": 0},
#         'Very smooth': {"D": 8.23, "E": 3.33, "F": -0.036},
#     }
#
#     exterior_heat_transfer_coefficient = material[material_roughness]["D"] + material[material_roughness]["E"] * local_wind_speed + material[material_roughness]["F"] * np.power(local_wind_speed, 2)
#     return exterior_heat_transfer_coefficient
#
#
# def generate_numerous_vectors(samples=1000):
#     """
#     Create a set of rays from a point source
#
#     Parameters
#     ----------
#     samples : int
#         Number of rays to generate
#
#     Returns
#     -------
#     vector : array(x, y, z)
#         Vector direction for each ray
#     theta : array(float)
#         Angle (in radians) between point source plane and vector end-point
#     sky : array(bool)
#         denotes whether ray cast is above/below ground (True is above ground - towards the sky)
#     """
#
#     offset = 2 / samples
#     increment = np.pi * (3 - np.sqrt(5))
#     y = ((np.arange(samples) * offset) - 1) + (offset / 2)
#     r = np.sqrt(1 - np.power(y, 2))
#     phi = np.arange(samples) * increment
#     x = np.cos(phi) * r
#     z = np.sin(phi) * r
#     vector = np.array([x, y, z]).T
#     theta = np.fabs(np.arctan(z / np.sqrt(np.power(x, 2) + np.power(y, 2))))
#     sky = z > 0
#     return vector, theta, sky
#
#
# def sky_view_factor(altitude, is_sky):
#     """
#     Calculate the sky view factor between source and sample vector (including below horizon effects.
#
#     Source: Jianxiang Huang, Jose Guillermo Cedeño-Laurent, John D. Spengler. 2014. CityComfort+: A simulation-based method for predicting mean radiant temperature in dense urban areas, https://doi.org/10.1016/j.buildenv.2014.05.019
#
#     Parameters
#     ----------
#     altitude : float
#         Angle between horizon and sample point (radians)
#     is_sky : bool
#         Different method used for below horizontal samples (accounting for person height)
#
#     Returns
#     -------
#     sky_view_factor : array(float)
#         Sky view factor value
#     """
#     az = np.pi / 4  # 45 degrees
#     sky_view_factor = np.where(is_sky, 0.0355 * np.sin(altitude) + 2.33 * np.cos(altitude) * np.sqrt(
#         0.0213 * np.power(np.cos(az), 2) + 0.0091 * np.power(np.sin(az), 2)), 0)
#     return sky_view_factor
#
#
# def radiation_on_square_metre(sample_vectors_altitude, sample_vectors_is_sky, sample_vector_radiation):
#     A = sample_vectors_altitude[sample_vectors_is_sky]
#     B = np.array([i[sample_vectors_is_sky] for i in sample_vector_radiation])
#     C = np.sin(A) * B
#     radiation = C.sum(axis=1)
#     return radiation
#
#
# def surface_temperature(dry_bulb_temperature, ground_temperature, ground_1m2_radiation, ground_emissivity,
#                         ground_absorptivity, ground_k_value, ground_thickness,
#                         ground_convective_heat_transfer_coefficient):
#     a = ground_emissivity * STEFAN_BOLTZMANN_CONSTANT
#     b = ground_k_value / ground_thickness + ground_convective_heat_transfer_coefficient
#     c = -(ground_k_value * (
#             ground_temperature + KELVIN) / ground_thickness + ground_1m2_radiation * ground_absorptivity + ground_convective_heat_transfer_coefficient * (
#                   dry_bulb_temperature + KELVIN))
#
#     Ts = []
#     X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
#             np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                               np.sqrt(3) * np.sqrt(
#                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                               2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
#         (4 * (2 / 3) ** (1 / 3) * c) / (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                 2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
#             2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
#     if X.sum() > 0:
#         Ts.append(X)
#
#     X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
#             np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                               np.sqrt(3) * np.sqrt(
#                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                               2 ** (1 / 3) * 3 ** (2 / 3) * a)) + 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
#         (4 * (2 / 3) ** (1 / 3) * c) / (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                 2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
#             2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
#     if X.sum() > 0:
#         Ts.append(X)
#
#     X = -0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
#             np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                                np.sqrt(3) * np.sqrt(
#                            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                                2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt((2 * b) / (a * np.sqrt(
#         (4 * (2 / 3) ** (1 / 3) * c) / (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                 np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                 2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
#             2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
#         27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
#     if X.sum() > 0:
#         Ts.append(X)
#
#     X = 0.5 * np.sqrt((2 * b) / (a * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
#             np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                                                      np.sqrt(3) * np.sqrt(
#                                                  27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
#                                                      1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (
#                               np.sqrt(3) * np.sqrt(
#                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
#                               2 ** (1 / 3) * 3 ** (2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (
#                               np.sqrt(3) * np.sqrt(
#                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
#                               1 / 3)) - 1 / 2 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
#             np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
#                                                                 np.sqrt(3) * np.sqrt(
#                                                             27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
#                                                                 1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))
#     if X.sum() > 0:
#         Ts.append(X)
#
#     return Ts[0]
#
#
#
