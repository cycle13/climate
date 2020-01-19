import pathlib
import platform
import subprocess
import pandas as pd
import numpy as np
from pvlib.solarposition import get_solarposition
from scipy import spatial

from .constants import REINHART_PATCH_CONVERSION_FACTOR, TREGENZA_PATCH_CONVERSION_FACTOR
from .helpers import chunk


class SkyMatrix(object):

    def __init__(self, m_type="Reinhart"):
        self.patch_vectors = None
        self._type = None

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

def sun_position(datetime, latitude, longitude):
    """
    Calculate the solar position based on latitude, longitude and time-zone aware datetime.

    Parameters
    ----------
    datetime
        A datetime or list of datetimes
    latitude
        The latitude of the location being assessed
    longitude
        The longitude of the location being assessed

    Returns
    -------
    DataFrame
        Contains solar position information for the datetimes passed

    """
    solar_metrics = get_solarposition(datetime, latitude, longitude)
    solar_metrics.rename(columns={'apparent_zenith': 'solar_apparent_zenith_angle', 'zenith': 'solar_zenith_angle',
                                  'apparent_elevation': 'solar_apparent_elevation_angle',
                                  'elevation': 'solar_elevation_angle', 'azimuth': 'solar_azimuth_angle',
                                  'equation_of_time': 'solar_equation_of_time'}, inplace=True)
    return solar_metrics


def gendaymtx(wea_file, direct=True, reinhart=False):
    """
    Run the Radiance Gendaymtx program from an input WEA file.

    *These calculations are dependant on the excellent **Radiance gendaymtx** program, which is expected to be
    located in either "C:/Radiance/bin/gendaymtx" or "/usr/local/radiance/bin/gendaymtx"*

    Parameters
    ----------
    wea_file : basestring
        Path to WEA file containing direct and diffuse radiation data
    direct : bool
        True: Direct-radiation only, False: Diffuse-radiation only
    reinhart : bool
        True: Create a sky-matrix with 577 patches, False: Create a sky-matrix with 146 patches

    Returns
    -------
    radiation_matrix : np.array
        A Numpy matrix of shape 8760*n, corresponding to each hour of the year, and each patch value

    """

    # Create output file path
    if direct:
        mtx_file = pathlib.Path(wea_file).with_suffix(".dirmtx")
    else:
        mtx_file = pathlib.Path(wea_file).with_suffix(".diffmtx")

    # Create run command
    if platform.system() != "Windows":
        cmd = '"/usr/local/radiance/bin/gendaymtx" -m {0:} -{1:} -O1 "{2:}" > "{3:}"'.format(2 if reinhart else 1,
                                                                                             "d" if direct else "s",
                                                                                             wea_file,
                                                                                             mtx_file)
    else:
        cmd = '"C:/Radiance/bin/gendaymtx" -m {0:} -{1:} -O1 "{2:}" > "{3:}"'.format(2 if reinhart else 1,
                                                                                     "d" if direct else "s",
                                                                                     wea_file, mtx_file)

    # Run command
    subprocess.run(cmd, shell=True)

    # Load the resultant annual patch-value matrix
    radiation_matrix = load_sky_matrix(mtx_file, reinhart)

    return radiation_matrix


def load_sky_matrix(mtx_file, reinhart=False):
    """
    Load an existing matrix file created by gendaymtx.

    Parameters
    ----------
    mtx_file : string
        Path to matrix file containing either direct or diffuse radiation
    reinhart : bool
        True: Load a sky-matrix with 577 patches, False: Load a sky-matrix with 146 patches

    Returns
    -------
    radiation_matrix : np.array
        A Numpy matrix of shape 8760*n, corresponding to each hour of the year, and each patch value

    """
    radiation_matrix = pd.read_csv(mtx_file, sep="\s+", skip_blank_lines=True, skiprows=8, header=None).values
    radiation_matrix = np.sum(radiation_matrix * np.array([0.265074126, 0.670114631, 0.064811243]), axis=1)
    radiation_matrix = np.array(list(chunk(radiation_matrix, n=8760, method="size"))[1:])
    radiation_matrix = np.multiply(radiation_matrix.T, REINHART_PATCH_CONVERSION_FACTOR if reinhart else TREGENZA_PATCH_CONVERSION_FACTOR)
    return radiation_matrix


def create_sky_matrices(wea_file, reinhart=False):
    """
    Create sky-matrices for the given WEA file.

    *These calculations are dependant on the excellent **Radiance gendaymtx** program, which is expected to be
    located in either "C:/Radiance/bin/gendaymtx" or "/usr/local/radiance/bin/gendaymtx"*

    Parameters
    ----------
    wea_file : basestring
        Path to WEA file containing direct and diffuse radiation data
    direct : bool
        True: Direct-radiation only, False: Diffuse-radiation only
    reinhart : bool
        True: Create a sky-matrix with 577 patches, False: Create a sky-matrix with 146 patches

    Returns
    -------
    sky_matrix : np.array
        A Numpy matrix of shape 8760*n, corresponding to each hour of the year, and each patch value

    """

    direct_sky_matrix = gendaymtx(wea_file, direct=True, reinhart=reinhart)
    print("Direct sky matrix calculated: {0:}".format(pathlib.Path(wea_file).with_suffix(".dirmtx")))

    diffuse_sky_matrix = gendaymtx(wea_file, direct=False, reinhart=reinhart)
    print("Diffuse sky matrix calculated: {0:}".format(pathlib.Path(wea_file).with_suffix(".diffmtx")))

    total_sky_matrix = direct_sky_matrix + diffuse_sky_matrix

    return direct_sky_matrix, diffuse_sky_matrix, total_sky_matrix
