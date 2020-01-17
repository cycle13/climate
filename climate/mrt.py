import numpy as np
from .environment import sky_emissivity
from .material import OpaqueMaterial

# from .constants import KELVIN, STEFAN_BOLTZMANN_CONSTANT

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
    return vector


def vector_horizon_angle(vector):
    """
    Returns the angle between a given vector or vectors, and the horizon. +Ve for vectors above horizon

    Parameters
    ----------
    vector : ndarray
        3 dimensional vector or array of vectors

    Returns
    -------
    angle : float
        Angle between input vector/s and z0 plane
    """
    vector = np.array(vector)
    if len(vector.shape) == 1:
        return np.arctan(vector[2] / np.sqrt(np.power(vector[0], 2) + np.power(vector[1], 2)))
    else:
        angles = []
        for vec in vector:
            angles.append(np.arctan(vec[2] / np.sqrt(np.power(vec[0], 2) + np.power(vec[1], 2))))
        return np.array(angles)


def sky_view_factor(vectors, shading_geometry=None):
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

    if shading_geometry is not None:
        raise NotImplementedError("Not yet implemented!")

    az = np.pi / 4  # 45 degrees
    altitude = vector_horizon_angle(vectors)
    sky_view_factor = np.where(altitude > 0, 0.0355 * np.sin(altitude) + 2.33 * np.cos(altitude) * np.sqrt(
        0.0213 * np.power(np.cos(az), 2) + 0.0091 * np.power(np.sin(az), 2)), 0)
    return sky_view_factor


def incident_radiation_at_angle(radiation, angle, degrees=False):
    return np.sin(np.radians(angle)) * radiation if degrees else np.sin(angle) * radiation

# TODO - Add method for squuare metre radiation calculating incidernt radiaiton total from all sky patch vectors and readaiton values


def radiation_on_square_metre(sample_vectors_altitude, sample_vectors_is_sky, sample_vector_radiation):
    A = sample_vectors_altitude[sample_vectors_is_sky]
    B = np.array([i[sample_vectors_is_sky] for i in sample_vector_radiation])
    C = np.sin(A) * B
    radiation = C.sum(axis=1)
    return radiation
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
