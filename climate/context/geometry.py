import uuid
import numpy as np
from scipy.spatial import Delaunay

from climate.context.material import Material


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
