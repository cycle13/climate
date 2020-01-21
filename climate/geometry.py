import numpy as np
from scipy.spatial import KDTree
from scipy.interpolate import LinearNDInterpolator


# Point/vector methods


def fibonacci_sphere(n_points=1000, cartesian=True):
    """ Create a set of evenly distributed points around a sphere

    Parameters
    ----------
    n_points : int
        Number of points to generate

    Returns
    -------
    points : nDArray
        Location of sample point on sphere
    """
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
    if cartesian:
        return np.stack([x, y, z]).T
    else:
        return np.stack([theta, phi, [1] * n_points])


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


def closest_point(sample_points, existing_points, n_closest=1):
    """ Find the n-closest points within a set of existing points

    Parameters
    ----------
    sample_points : array(x, y, z)
        Points to sample for
    existing_points : array(x, y, z)
        Set of points which will be sampled from
    n_closest : int
        The number of "near" points to return

    Returns
    -------
    distances : array(float)
        Distance between each source point and nearest target point/s
    indices : array(int)
        Indices of each target point near to the source point/s

    """
    return KDTree(existing_points).query(sample_points, n_closest)


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


def resample_point_values(sample_points, existing_points, existing_values):
    """ Triangulate values at sample points from known point-value locations

    Parameters
    ----------
    sample_points : ndArray
        Point about which values will be interpolated
    existing_points : ndArray
        Known points
    existing_values : ndArray
        Known values

    Returns
    -------
    sample_values : ndArray
        Re-sampled (triangulated) point values
    """
    # TODO - Add checks to ensure existing points and values shapes match
    return LinearNDInterpolator(existing_points, existing_values)(sample_points).T


def view_factor(vector, shading_geometry=None):
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


