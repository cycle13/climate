import numpy as np
from scipy.spatial import KDTree, Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
import uuid
import pathlib
from .material import Material

class Point(object):
    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        self.guid = str(uuid.uuid4())
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Point({}, {}, {})".format(self.x, self.y, self.z)

    def __mul__(self, other: float):
        return Point(self.x * other, self.y * other, self.z * other)

    def unitize(self):
        """ Gets the unitized form of the point/vector

        Returns
        -------
        unitized_vector : Point
            Unitized form of point/vector
        """
        magnitude = np.sqrt(np.power(self.x, 2) + np.power(self.y, 2) + np.power(self.z, 2))
        return Point((self.x / magnitude, self.y / magnitude, self.z / magnitude))

    def array(self):
        return np.array([self.x, self.y, self.z])

class PointCollection(object):
    def __init__(self, points: np.ndarray):
        self.points = points
        self.points_x = [i.x for i in points]
        self.points_y = [i.y for i in points]
        self.points_z = [i.z for i in points]

    def array(self):
        return

class Polygon(object):
    def __init__(self, outer_vertices: np.ndarray = None, inner_vertices: np.ndarray = None, material: Material = Material(), inner_material: Material = None):
        self.guid = str(uuid.uuid4())
        self.outer_vertices = outer_vertices
        self.inner_vertices = inner_vertices
        self.tri = self.triangulate()
        self.material = material

    def __repr__(self):
        return "Polygon:\n- Material: {0:}\n- Tri: {1:} triangulated sub-polygons".format(self.material.guid, len(self.tri))

    def triangulate(self, inner: bool = False):
        return triangulate_3d_surfaces(self.outer_vertices, self.inner_vertices)

    def to_rad_string(self):
        rad_string = []
        rad_string.append(self.material.to_rad_string())
        for bp in self.tri:
            rad_string.append(rad_string_polygon(bp, id=self.guid, material=self.material.guid))
        return "\n\n".join(rad_string)


# Point/vector methods
def rad_string_polygon(boundary_points: np.ndarray, id: str=str(uuid.uuid4()), material: str="material_id"):
    return "{0:} polygon {1:}\n0\n0\n{2:} ".format(id, material, len(boundary_points) / 3) + " ".join([str(i) for i in boundary_points.flatten()])


def fibonacci_sphere(n_points: int=1000, cartesian: bool=True):
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


def vector_horizon_angle(vector: np.ndarray):
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


def closest_point(sample_points: np.ndarray, existing_points: np.ndarray, n_closest: int=1):
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


def unit_vector(start: np.ndarray, end: np.ndarray):
    """
    Returns the unit vector of a line described by its start and end points
    :type start: [x, y] coordinate array
    :type end: [x, y] coordinate array
    :return: [x, y] vector array
    """
    pt_distance = np.array(end) - np.array(start)
    vector = pt_distance / np.sqrt(np.sum(pt_distance * pt_distance))
    return vector


def normalise_vector(vector: np.ndarray):
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


def angle_between(v1: np.ndarray, v2: np.ndarray, degrees: bool=False):
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


def view_factor(vector: np.ndarray, shading_geometry=None):
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

def translated_point(uv: np.ndarray, uw: np.ndarray, origin: np.ndarray, point: np.ndarray):
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


def untranslated_point(uv: np.ndarray, uw: np.ndarray, origin: np.ndarray, point: np.ndarray):
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


def triangulate_3d_surfaces(parent_surface_vertices: np.ndarray, child_surfaces_vertices: np.ndarray=None):
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

    if child_surfaces_vertices is not None:
        child_surfaces_vertices_translated = np.array(
            [[translated_point(uv, uw, parent_surface_vertices[0], i) for i in ch] for ch in child_surfaces_vertices])
    else:
        child_surfaces_vertices_translated = None

    parent_points = parent_surface_vertices_translated
    child_points = [item for sublist in child_surfaces_vertices_translated for item in
                    sublist] if child_surfaces_vertices is not None else None

    points = np.concatenate([parent_points, child_points]) if child_surfaces_vertices is not None else parent_points
    tri = Delaunay(points).simplices.copy()

    if child_surfaces_vertices is not None:
        mask = []
        for face_pts in points[tri]:
            n = []
            for child_pts in child_surfaces_vertices_translated:
                n.append(len(np.array([x for x in set(tuple(x) for x in face_pts) & set(tuple(x) for x in child_pts)])))
            if 3 in n:
                mask.append(False)
            else:
                mask.append(True)
    else:
        mask = None

    triangulated_surface_vertices = []
    lm = points[tri][mask] if child_surfaces_vertices is not None else points[tri]
    for i in lm:
        mm = []
        for j in i:
            mm.append(untranslated_point(uv, uw, parent_surface_vertices[0], j))
        triangulated_surface_vertices.append(mm)

    return np.array(triangulated_surface_vertices)

# Geometry to Radiance methods
# Geometry to Radiance with Material methods
# Proper context object polygon thing

class Point(object):
    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        self.guid = str(uuid.uuid4())
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Point({}, {}, {})".format(self.x, self.y, self.z)

    def __mul__(self, other: float):
        return Point(self.x * other, self.y * other, self.z * other)

    def unitize(self):
        """ Gets the unitized form of the point/vector

        Returns
        -------
        unitized_vector : Point
            Unitized form of point/vector
        """
        magnitude = np.sqrt(np.power(self.x, 2) + np.power(self.y, 2) + np.power(self.z, 2))
        return Point(self.x / magnitude if magnitude > 0 else 0, self.y / magnitude if magnitude > 0 else 0,
                     self.z / magnitude if magnitude > 0 else 0)

    def array(self):
        return np.array([self.x, self.y, self.z])


class Points(object):
    def __init__(self, points: np.ndarray):
        self.points = np.array(points)

        self.xs = np.array([i.x for i in points])
        self.ys = np.array([i.y for i in points])
        self.zs = np.array([i.z for i in points])

    def xmin(self):
        return np.min(self.xs)

    def xmax(self):
        return np.max(self.xs)

    def ymin(self):
        return np.min(self.ys)

    def ymax(self):
        return np.max(self.ys)

    def zmin(self):
        return np.min(self.zs)

    def zmax(self):
        return np.max(self.zs)

    def xmedian(self):
        return np.median(self.xs)

    def ymedian(self):
        return np.median(self.ys)

    def zmedian(self):
        return np.median(self.zs)

    def centroid(self):
        return Point(self.xmedian(), self.ymedian(), self.zmedian())

    def __repr__(self):
        return "Points: {} points".format(len(self.points))

    def array(self):
        return np.vstack([self.xs, self.ys, self.zs]).T

    def unitize(self):
        return np.array([pt.unitize() for pt in self.points])


class SphericalCoordinate(object):
    def __init__(self, theta: float = 0, phi: float = 0):
        self.theta = theta
        self.phi = phi

    def __repr__(self):
        return "SphericalCoordinate: ({}, {})".format(self.theta, self.phi)

    def cartesian(self):
        x = np.cos(self.theta) * np.sin(self.phi)
        y = np.sin(self.theta) * np.sin(self.phi)
        z = np.cos(self.phi)
        return Point(x, y, z)


class SphericalCoordinates(object):
    def __init__(self, spherical_coordinates: np.ndarray):
        self.spherical_coordinates = np.array(spherical_coordinates)

        self.thetas = np.array([i.theta for i in spherical_coordinates])
        self.phis = np.array([i.phi for i in spherical_coordinates])

    def __repr__(self):
        return "SphericalCoordinates: {} spherical coordinates".format(len(self.spherical_coordinates))


#     def cartesian(self):
#         return np.array([sph.cartesian() for sph in self.spherical_coordinates])


def SkyPatchLocations(m=1):
    # Create the number of patches per row
    row_patches_base = np.array([1, 6, 12, 18, 24, 24, 30, 30])
    row_patches_count = np.concatenate([np.array([1]), np.repeat(row_patches_base[1:] * m, m)])

    # Generate azimuth angles for each patch
    row_azimuths = [[0]] + [np.linspace(0, 2 * np.pi, i, endpoint=False).tolist() for i in row_patches_count[1:]]

    # Generate altitude angles for each patch
    row_altitudes = np.linspace(0, np.pi / 2, len(row_patches_count)).tolist()

    # Construct the theta/phi coordinates array
    theta = np.hstack(row_azimuths)
    phi = np.hstack([[i] * row_patches_count[n] for n, i in enumerate(row_altitudes)])

    return SphericalCoordinates(theta, phi)

# Generate point locations on sphere

def sphere_cartesian(theta, phi):
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    return np.vstack([x, y, z]).T


def SkyPatchLocations(m=1):
    # Create the number of patches per row
    row_patches_base = np.array([1, 6, 12, 18, 24, 24, 30, 30])
    row_patches_count = np.concatenate([np.array([1]), np.repeat(row_patches_base[1:] * m, m)])

    # Generate azimuth angles for each patch
    row_azimuths = [[0]] + [np.linspace(0, 2 * np.pi, i, endpoint=False).tolist() for i in row_patches_count[1:]]

    # Generate altitude angles for each patch
    row_altitudes = np.linspace(0, np.pi / 2, len(row_patches_count)).tolist()

    # Construct the theta/phi coordinates array
    theta = np.hstack(row_azimuths)
    phi = np.hstack([[i] * row_patches_count[n] for n, i in enumerate(row_altitudes)])

    return sphere_cartesian(theta, phi)


# Construct sphere points using the spherical coordinates

#
# x, y, z = sphere_cartesian(theta, phi).T
# fig = plt.figure()
# ax = Axes3D(fig)
#
# ax.scatter(x, y, z)
# print(x.shape)
