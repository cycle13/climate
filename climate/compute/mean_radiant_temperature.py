from scipy.interpolate import bisplev

from climate.common.helpers import *
from climate.common.constants import ray_last_bounce_vector, ray_intersected_points

import numpy as np

# TODO - Add a rolling thermal lag to the calculated surface temperature for the ground - see added file for surface_tmeperature


def fanger_solar_exposure(solar_azimuth_angle, solar_elevation_angle, posture=1):
    # Sitting bivariate B-spline and its derivatives
    sitting_bisplrep = [[0, 0, 0, 0, 0, 180, 180, 180, 180, 180], [0, 0, 0, 0, 0, 90, 90, 90, 90, 90],
                        [0.42999558258469306, 0.49777802226651985, 0.2858541264803382, 0.2636331839991635,
                         0.10059901058304405, 0.5904998653021177, 0.6393605287969937, 0.41177803047742195,
                         0.16397939147762605, 0.1145630272512949, -0.07290688451711066, -0.0877360565501316,
                         0.03719879969518147, 0.06771059788029093, 0.09444998526069391, 0.5573351684449549,
                         0.6212235986152396, 0.3384990152297299, 0.25505266892999545, 0.1011441730110879,
                         0.4432956571996788, 0.49809124858382825, 0.29471168936411446, 0.19682482035937438,
                         0.10008130856803796], 4, 4]

    # Standing bivariate B-spline and its derivatives
    standing_bisplrep = [[0, 0, 0, 0, 0, 180, 180, 180, 180, 180], [0, 0, 0, 0, 0, 90, 90, 90, 90, 90],
                         [0.365433469329803, 0.41471995039390336, 0.3539202584010255, 0.35205668670475776,
                          0.21505967838173534, 0.5304745700779437, 0.6180584137541132, 0.11434859278048302,
                          0.4862162611010728, 0.20252438358996272, -0.015147290187610778, 0.22189948439503024,
                          0.6990946114268216, -0.000718703369787728, 0.22472889635480628, 0.5176922764465676,
                          0.35055123160310636, -0.0032935618498728487, 0.3404006983313149, 0.19936403473400507,
                          0.37870178660536147, 0.24613731159172733, 0.06300314787643235, 0.23364607863218287,
                          0.2171651821703637], 4, 4]

    # Convert inputs to 0-90, 0-180 range arrays
    azimuth = np.where(solar_azimuth_angle < 180, solar_azimuth_angle, 360 - solar_azimuth_angle)
    altitude = np.where(solar_elevation_angle < 0, 0, np.where(solar_elevation_angle > 90, 90 - solar_elevation_angle, solar_elevation_angle))

    # Check for equal length arrays
    if not len(azimuth) == len(altitude):
        raise ValueError("Inputs should all be the same length")

    if posture == 0:
        solar_exposure = np.array([bisplev(azimuth[i], altitude[i], sitting_bisplrep) for i in range(len(azimuth))])
    elif posture == 1:
        solar_exposure = np.array([bisplev(azimuth[i], altitude[i], standing_bisplrep) for i in range(len(azimuth))])
    else:
        solar_exposure = None

    return solar_exposure


def solar_adjusted_mean_radiant_temperature(dry_bulb_temperature, diffuse_horizontal_radiation, global_horizontal_radiation, direct_normal_radiation, solar_azimuth, solar_altitude, clothing_absorbtivity=0.7, ground_reflectivity=0.4, shading_transmissivity=1.0):
    fraction_body_exposed_radiation = 0.71  # Fraction of the body visible to radiation. 0.725 is for standing, 0.696 for sitting, 0.68 for lying down. A nominal value of 0.71 is used.
    radiative_heat_transfer_coefficient = 6.012  # A good guess at the radiative heat transfer coefficient
    projected_area_fraction = fanger_solar_exposure(solar_azimuth, solar_altitude, posture=1)

    effective_radiant_flux = ((0.5 * fraction_body_exposed_radiation * (
            diffuse_horizontal_radiation + (global_horizontal_radiation * ground_reflectivity)) + (
                                       fraction_body_exposed_radiation * projected_area_fraction * direct_normal_radiation)) * shading_transmissivity) * (
                                     clothing_absorbtivity / 0.95)

    mean_radiant_temperature_delta = (
            effective_radiant_flux / (fraction_body_exposed_radiation * radiative_heat_transfer_coefficient))

    solar_adjusted_mrt = dry_bulb_temperature + mean_radiant_temperature_delta

    return solar_adjusted_mrt


def mrt_solar_adjusted(self, clothing_absorbtivity=0.7, ground_reflectivity=0.4, shading_transmissivity=1.0):

    mrt_sa = solar_adjusted_mean_radiant_temperature(
        self.dry_bulb_temperature, self.diffuse_horizontal_radiation, self.global_horizontal_radiation, self.direct_normal_radiation, self.solar_azimuth_angle,
        self.solar_elevation_angle, clothing_absorbtivity=clothing_absorbtivity, ground_reflectivity=ground_reflectivity,
        shading_transmissivity=shading_transmissivity)
    self.mean_radiant_temperature_solar_adjusted = pd.Series(index=self.index, name="mean_radiant_temperature_solar_adjusted", data=mrt_sa)
    print("Mean radiant temperature (solar adjusted) calculations successful")
    return self.mean_radiant_temperature_solar_adjusted


def calc_sky_emissivity(dew_point_temperature, total_sky_cover):
    """
    Calculate the sky emissivity using the Clark-Allen method used by EnergyPlus.

    Clark and Allen report that the standard estimation error for the atmospheric long-wave radiation was 10 𝑊𝑊/𝑚𝑚2.

    Parameters
    ----------
    dew_point_temperature : float
        Dew-point temperature of air in degrees C
    total_sky_cover : int
        Total sky cover in tenths. Conversion to decimal 0-1 happens within the method.

    Returns
    -------
    sky_emissivity : float

    """

    sky_emissivity = (0.787 + 0.764 * np.log((dew_point_temperature + KELVIN) / KELVIN)) * (
            1 + 0.0224 * total_sky_cover - 0.0035 * np.power(total_sky_cover, 2) + 0.0028 * np.power(
        total_sky_cover, 3))
    return sky_emissivity


def simple_combined_exterior_heat_transfer_coefficient(material_roughness, local_wind_speed):
    """
    The simple algorithm uses surface roughness and local surface windspeed to calculate the exterior heat transfer coefficient.

    Note that the simple correlation yields a combined convection and radiation heat transfer coefficient.

    Material roughness correlations are taken from Figure 1, Page 22.4, ASHRAE Handbook of Fundamentals (ASHRAE 1989).

    Source: http://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-020.html for further details.

    Parameters
    ----------
    material_roughness : string
        A surface roughness index. Acceptable values are: "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth" and "Very smooth"
    local_wind_speed : float
        Local wind speed calculated at the height above ground of the surface centroid

    Returns
    -------
    heat_transfer_coefficient
        The surface heat transfer coefficient for the variables passed
    """

    material = {
        "Very rough": {"D": 11.58, "E": 5.89, "F": 0},
        "Rough": {"D": 12.49, "E": 4.065, "F": 0.028},
        'Medium rough': {"D": 10.79, "E": 4.192, "F": 0},
        'Medium smooth': {"D": 8.23, "E": 4, "F": -0.057},
        'Smooth': {"D": 10.22, "E": 3.1, "F": 0},
        'Very smooth': {"D": 8.23, "E": 3.33, "F": -0.036},
    }

    exterior_heat_transfer_coefficient = material[material_roughness]["D"] + material[material_roughness]["E"] * local_wind_speed + material[material_roughness]["F"] * np.power(local_wind_speed, 2)
    return exterior_heat_transfer_coefficient


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
    theta = np.fabs(np.arctan(z / np.sqrt(np.power(x, 2) + np.power(y, 2))))
    sky = z > 0
    return vector, theta, sky


def sky_view_factor(altitude, is_sky):
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
    az = np.pi / 4  # 45 degrees
    sky_view_factor = np.where(is_sky, 0.0355 * np.sin(altitude) + 2.33 * np.cos(altitude) * np.sqrt(
        0.0213 * np.power(np.cos(az), 2) + 0.0091 * np.power(np.sin(az), 2)), 0)
    return sky_view_factor


def radiation_on_square_metre(sample_vectors_altitude, sample_vectors_is_sky, sample_vector_radiation):
    A = sample_vectors_altitude[sample_vectors_is_sky]
    B = np.array([i[sample_vectors_is_sky] for i in sample_vector_radiation])
    C = np.sin(A) * B
    radiation = C.sum(axis=1)
    return radiation


def surface_temperature(dry_bulb_temperature, ground_temperature, ground_1m2_radiation, ground_emissivity,
                        ground_absorptivity, ground_k_value, ground_thickness,
                        ground_convective_heat_transfer_coefficient):
    a = ground_emissivity * STEFAN_BOLTZMANN_CONSTANT
    b = ground_k_value / ground_thickness + ground_convective_heat_transfer_coefficient
    c = -(ground_k_value * (
            ground_temperature + KELVIN) / ground_thickness + ground_1m2_radiation * ground_absorptivity + ground_convective_heat_transfer_coefficient * (
                  dry_bulb_temperature + KELVIN))

    Ts = []
    X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a)) + 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = -0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                               np.sqrt(3) * np.sqrt(
                           27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                               2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt((2 * b) / (a * np.sqrt(
        (4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
            2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
        27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
    if X.sum() > 0:
        Ts.append(X)

    X = 0.5 * np.sqrt((2 * b) / (a * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                     np.sqrt(3) * np.sqrt(
                                                 27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                     1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                              2 ** (1 / 3) * 3 ** (2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (
                              np.sqrt(3) * np.sqrt(
                          27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                              1 / 3)) - 1 / 2 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
            np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                                np.sqrt(3) * np.sqrt(
                                                            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                                1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))
    if X.sum() > 0:
        Ts.append(X)

    return Ts[0]


def openfield_mean_radiant_temperature(self, ground_roughness="Medium rough", ground_emissivity=0.8, ground_absorptivity=0.6, shading_transmissivity=1):

    # Sense check
    if ground_emissivity <= 0 or ground_emissivity > 1:
        raise Exception('ground_emissivity should be > 0 and < 1')

    if ground_roughness not in ["Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth", "Very smooth"]:
        raise Exception('ground_roughness should be one of "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth", "Very smooth"')

    if shading_transmissivity <0 or shading_transmissivity > 1:
        raise Exception('shading_transmissivity should be >= 0 and <1')

    if ground_absorptivity <= 0 or ground_absorptivity > 1:
        raise Exception('ground_absorptivity should be > 0 and < 1')

    ground_albedo = 1 - ground_emissivity
    ground_thickness = 1
    sky_emissivity = calc_sky_emissivity(self.dew_point_temperature, self.total_sky_cover / 10)
    ground_convective_heat_transfer_coefficient = simple_combined_exterior_heat_transfer_coefficient(ground_roughness, self.pedestrian_wind_speed)

    ground_k_value = np.interp(np.power(self.relative_humidity, 3), [0, 1e6], [0.33, 1.4])  # TODO: Check what this value is and where it comes from!

    ground_temperature = self.ground_temperature_500_weatherfile

    # Generate sample vectors for radiation from the sky
    sample_vectors, sample_vectors_altitude, sample_vectors_is_sky = generate_numerous_vectors(samples=1000)
    self.sample_vectors = sample_vectors
    self.sample_vectors_altitude = sample_vectors_altitude
    self.sample_vectors_is_sky = sample_vectors_is_sky

    # Find the closest patch value indices and distances from the sample vectors
    sample_vector_closest_patch_vector_distances, sample_vector_closest_patch_vector_indices = closest_point(sample_vectors, self.patch_centroids, n_closest=3)

    # Calculate the sky view factors for each of the sample vectors.
    sky_view_factors = sky_view_factor(sample_vectors_altitude, sample_vectors_is_sky) * shading_transmissivity  # Shade adjusted MRT effects (this modifies the visibility of the sky for both radiation received, and lost to sky)

    # Calculate the closest points in the sky matrix to the rays sent from the sample point
    _, ray_last_bounce_vector_closest_sample_vector_indices = closest_point(ray_last_bounce_vector, sample_vectors, n_closest=1)

    # Calculate the sky radiation values received by each sample vector
    sample_vector_radiation = []
    for total_sky_matrix_hour in self.total_sky_matrix:
        n_values = (total_sky_matrix_hour[sample_vector_closest_patch_vector_indices] * sample_vector_closest_patch_vector_distances).sum(axis=1) / sample_vector_closest_patch_vector_distances.sum(axis=1)
        n_values = np.where(sample_vectors[:, 2] <= 0, 0, n_values)  # Replace values where vectors are below ground
        total_sky_matrix_hour_radiation = sum(total_sky_matrix_hour)  # Total hourly radiation from original sky matrix
        resampled_sky_matrix_hour_radiation = sum(n_values)  # Total hourly radiation from sampled sky matrix
        sample_vector_radiation.append(np.where(n_values != 0, n_values / resampled_sky_matrix_hour_radiation * total_sky_matrix_hour_radiation, 0))
    sample_vector_radiation = np.array(sample_vector_radiation)

    # Calculate the total radiation received on 1 square metre
    ground_1m2_radiation = radiation_on_square_metre(sample_vectors_altitude, sample_vectors_is_sky, sample_vector_radiation)

    # Calculate radiation coming from the sky matrix (incorporating the view factor)
    sky_sourced_radiation = (sky_view_factors[ray_last_bounce_vector_closest_sample_vector_indices] * sample_vector_radiation[:, ray_last_bounce_vector_closest_sample_vector_indices])

    # Calculate radiation reflected by the ground
    ground_reflected_sky_radiation = (np.power(ground_albedo, ray_intersected_points) * sky_sourced_radiation).sum(axis=1)

    # Calculate ground surface temperature
    ground_surface_temperature = surface_temperature(self.dry_bulb_temperature, ground_temperature, ground_1m2_radiation, ground_emissivity,
                                                     ground_absorptivity, ground_k_value, ground_thickness,
                                                     ground_convective_heat_transfer_coefficient)

    # Calculate atmospheric radiation
    SVF = 1
    atmospheric_sky_radiation = sky_emissivity * self.horizontal_infrared_radiation_intensity * SVF

    # Calculate solar radiation
    ap = 0.7
    atmospheric_solar_radiation = ground_reflected_sky_radiation * ap

    # Calculate ground reflected/infrared radiation
    A1 = 1.2
    A2 = 1
    ground_radiation = (A1 / A2) * ground_emissivity * (STEFAN_BOLTZMANN_CONSTANT * ground_emissivity * np.power(ground_surface_temperature, 4)) * 0.4

    # Calculate experienced Mean Radiant Temperature
    mean_radiant_temperature = np.power(((atmospheric_sky_radiation + atmospheric_solar_radiation + ground_radiation) / STEFAN_BOLTZMANN_CONSTANT), 0.25) - KELVIN
    return mean_radiant_temperature


def mrt_openfield(self, ground_roughness="Medium rough", ground_emissivity=0.8, ground_absorptivity=0.6, shading_transmissivity=1):
    mrt_of = openfield_mean_radiant_temperature(self, ground_roughness=ground_roughness, ground_emissivity=ground_emissivity, ground_absorptivity=ground_absorptivity, shading_transmissivity=shading_transmissivity)
    mrt_of.index = self.index
    mrt_of.name = "mean_radiant_temperature_openfield"
    self.mean_radiant_temperature_openfield = mrt_of
    print("Mean radiant temperature (openfield) calculations successful")
    return self.mean_radiant_temperature_openfield

