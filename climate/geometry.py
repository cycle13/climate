import numpy as np
from scipy.spatial import KDTree, Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

class Polygon(object):
    def __init__(self):
        self.vertices = None


class Geometry(object):
    def __init__(self, Material):
        self.yes = None

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


def unit_vector(start, end):
    """
    Returns the unit vector of a line described by its start and end points
    :type start: [x, y] coordinate array
    :type end: [x, y] coordinate array
    :return: [x, y] vector array
    """
    pt_distance = np.array(end) - np.array(start)
    vector = pt_distance / np.sqrt(np.sum(pt_distance * pt_distance))
    return vector


def normalise_vector(vector):
    """
    Returns the angle between vectors 'v1' and 'v2'

    Parameters
    ----------
    vector : vector
        Non-unitized n-dimensional vector

    Returns
    -------
    normalise_vector : ndarray
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
    v1_u = normalise_vector(v1)
    v2_u = normalise_vector(v2)
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

## TRIANGULATION

def translated_point(uv, uw, origin, point):
    """
    Translates a 3D point into a 2D point
    :type uv: array
    :type uw: array
    :type origin: array
    :type point: array
    :return: array
    """
    x = (point[0] - origin[0]) * uv[0] + (point[1] - origin[1]) * uv[1] + (point[2] - origin[2]) * uv[2]
    y = (point[0] - origin[0]) * uw[0] + (point[1] - origin[1]) * uw[1] + (point[2] - origin[2]) * uw[2]
    return x, y


def untranslated_point(uv, uw, origin, point):
    """
    Translates a 2D point into a 3D point
    :type uv: array
    :type uw: array
    :type origin: array
    :type point: array
    :return: array
    """
    x = origin[0] + uv[0] * point[0] + uw[0] * point[1]
    y = origin[1] + uv[1] * point[0] + uw[1] * point[1]
    z = origin[2] + uv[2] * point[0] + uw[2] * point[1]
    return x, y, z


def triangulate_3d_surfaces(parent_surface_vertices, child_surfaces_vertices):
    """
    Returns a set of vertices describing the delaunay mesh of a parent surface minus the child surfaces
    :type parent_surface_vertices: array
    :type child_surfaces_vertices: array
    :return: array
    """
    uv = unit_vector(parent_surface_vertices[0], parent_surface_vertices[1])
    uw = unit_vector(parent_surface_vertices[0], parent_surface_vertices[3])

    parent_surface_vertices_translated = np.array(
        [translated_point(uv, uw, parent_surface_vertices[0], i) for i in parent_surface_vertices])
    child_surfaces_vertices_translated = np.array(
        [[translated_point(uv, uw, parent_surface_vertices[0], i) for i in ch] for ch in child_surfaces_vertices])

    parent_points = parent_surface_vertices_translated
    child_points = [item for sublist in child_surfaces_vertices_translated for item in sublist]

    points = np.concatenate([parent_points, child_points])
    tri = Delaunay(points).simplices.copy()

    mask = []
    for face_pts in points[tri]:
        n = []
        for child_pts in child_surfaces_vertices_translated:
            n.append(len(np.array([x for x in set(tuple(x) for x in face_pts) & set(tuple(x) for x in child_pts)])))
        if 3 in n:
            mask.append(False)
        else:
            mask.append(True)

    triangulated_surface_vertices = []
    for i in points[tri][mask]:
        mm = []
        for j in i:
            mm.append(untranslated_point(uv, uw, parent_surface_vertices[0], j))
        triangulated_surface_vertices.append(mm)

    return np.array(triangulated_surface_vertices)

# Geometry to Radiance methods
# Geometry to Radiance with Material methods
# Proper geometry object polygon thing

