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


