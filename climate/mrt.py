import numpy as np
from scipy import spatial

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

    return spatial.KDTree(target_points).query(source_points, n_closest)


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


def sky_view_factor(vector, shading_geometry=None):
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
    altitude = vector_horizon_angle(vector)
    return np.where(altitude > 0, 0.0355 * np.sin(altitude) + 2.33 * np.cos(altitude) * np.sqrt(
        0.0213 * np.power(np.cos(az), 2) + 0.0091 * np.power(np.sin(az), 2)), 0)


def incident_radiation_at_angle(radiation, angle, degrees=False):
    """
    Returns radiation reduced by angle of incidence

    Parameters
    ----------
    radiation : float
        Radiation flux in W/m2
    angle : float
        Angle of incidence (in either degrees or radians in conjunction with flag)
    degrees : bool
        True if angle passed in in degrees, False if angle passed is in radians

    Returns
    -------
    Factored input radiance by angle of incidence reduction
    """
    return np.sin(np.radians(angle)) * radiation if degrees else np.sin(angle) * radiation

# TODO - Add methods
#  - Sky matrix resampling (based on new point layout - possibly add this to the SKY file instead as resampling method. Arguments include "old_patch_centroids",  "new_patch_centroids" and "n_samples"
#  - Total square metre radiation calculating incidernt radiaiton total from all sky patch vectors and radaiton values


# def resample_sky_matrix(sky_matrix, source_points, target_points, n_closest=1):
#     nv_dist, nv_idx = closest_point(source_points, target_points, n_closest=n_closest)
#     sample_vector_radiation = []
#     for total_sky_matrix_hour in sky_matrix:
#         n_values = (total_sky_matrix_hour[nv_idx] * nv_dist).sum(axis=1) / nv_dist.sum(axis=1)
#         n_values = np.where(target_points[:, 2] <= 0, 0, n_values)  # Replace values where vectors are below ground
#         total_sky_matrix_hour_radiation = sum(total_sky_matrix_hour)  # Total hourly radiation from original sky matrix
#         resampled_sky_matrix_hour_radiation = sum(n_values)  # Total hourly radiation from sampled sky matrix
#         sample_vector_radiation.append(
#             np.where(n_values != 0, n_values / resampled_sky_matrix_hour_radiation * total_sky_matrix_hour_radiation,
#                      0))
#     return np.array(sample_vector_radiation)

def radiation_flux():
    return 1


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
