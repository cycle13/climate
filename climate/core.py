from climate.constants import *
from climate.helpers import *
from climate.psychrometrics import psychrometric_calculations
from climate.sun import sun_position, sky_matrix_calculations, load_sky_matrix
from climate.environment import *
from climate.comfort import *

import platform
import subprocess

import pandas as pd
import numpy as np
import pathlib
from io import StringIO
# from climate.helpers import slugify, chunk, angle_between, renamer, unit_vector, wind_speed_at_height, ground_temperature_at_depth
from pvlib.solarposition import get_solarposition
from psychrolib import GetHumRatioFromRelHum
from scipy import spatial
from scipy.interpolate import bisplev

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from matplotlib import cm


class Weather(object):

    def __init__(self, file_path):
        """
        A weather object containing a representative years worth of climate data

        Parameters
        ----------
        file_path : string
            The path to the EPW file to be loaded and processed

        Returns
        -------
        weather
            A weather-object giving access to all loaded, calculated and assigned data

        """

        # Metadata
        self.file_path = pathlib.Path(file_path).absolute()

        # Location variables
        self.city = None
        self.region = None
        self.country = None
        self.dataset_type = None
        self.station_id = None
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.time_zone = None
        self.design_conditions = None
        self.typical_extreme_periods = None
        self.ground_temperatures = None
        self.holidays_daylight_savings = None
        self.comments_1 = None
        self.comments_2 = None

        # Series variables
        self.index = DATETIME_INDEX
        self.df = None
        self.data_source_and_uncertainty_flags = None
        self.dry_bulb_temperature = None
        self.dew_point_temperature = None
        self.relative_humidity = None
        self.atmospheric_station_pressure = None
        self.extraterrestrial_horizontal_radiation = None
        self.extraterrestrial_direct_normal_radiation = None
        self.horizontal_infrared_radiation_intensity = None
        self.global_horizontal_radiation = None
        self.direct_normal_radiation = None
        self.diffuse_horizontal_radiation = None
        self.global_horizontal_illuminance = None
        self.direct_normal_illuminance = None
        self.diffuse_horizontal_illuminance = None
        self.zenith_luminance = None
        self.wind_direction = None
        self.wind_speed = None
        self.total_sky_cover = None
        self.opaque_sky_cover = None
        self.visibility = None
        self.ceiling_height = None
        self.present_weather_observation = None
        self.present_weather_codes = None
        self.precipitable_water = None
        self.aerosol_optical_depth = None
        self.snow_depth = None
        self.days_since_last_snowfall = None
        self.albedo = None
        self.liquid_precipitation_depth = None
        self.liquid_precipitation_quantity = None
        self.ground_temperature_1 = None
        self.ground_temperature_2 = None
        self.ground_temperature_3 = None

        # Path variables
        self.wea_file = None
        self.csv_file = None

        # Calculated pedestrian wind speed
        self.pedestrian_wind_speed = None

        # Calculated sun-position variables
        self.solar_apparent_zenith_angle = None
        self.solar_zenith_angle = None
        self.solar_apparent_elevation_angle = None
        self.solar_elevation_angle = None
        self.solar_azimuth_angle = None
        self.solar_equation_of_time = None

        # Calculated psychrometric variables
        self.humidity_ratio = None
        self.wet_bulb_temperature = None
        self.partial_vapour_pressure_moist_air = None
        self.enthalpy = None
        self.specific_volume_moist_air = None
        self.degree_of_saturation = None

        # Calculated sky-dome variables
        self.direct_sky_matrix_path = None
        self.diffuse_sky_matrix_path = None
        self.reinhart = None
        self.direct_sky_matrix = None
        self.diffuse_sky_matrix = None
        self.total_sky_matrix = None
        self.patch_centroids = None
        self.patch_vectors = None
        self.patch_count = None

        # NV method values
        self.nv_sample_vectors = None
        self.mean_radiant_temperature_of = None
        self.mean_radiant_temperature_sa = None

        self.utci_sa = None
        self.utci_of = None

    # Methods below here for read/write
    def read(self, sky_matrix=False, reuse_matrix=True):
        """
        Read EPW weather-file into weather object.

        Additional processing is possible using the flags provided.

        Parameters
        ----------
        sky_matrix : bool
            Calculate the sky-dome, with each patch corresponding to an annual hourly set of direct, diffuse and
            total-radiation values.

            *These calculations are dependant on the excellent **Radiance gendaymtx** program, available from
            https://github.com/NREL/Radiance/releases/tag/5.2.*

        Returns
        -------
        weather
            A weather-object giving access to all the loaded and calculated data

        """
        with open(self.file_path, "r") as f:
            dat = f.readlines()

        # Read location data
        self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = dat[0].strip().split(",")[1:]
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)
        self.time_zone = float(self.time_zone)
        self.elevation = float(self.elevation)
        self.design_conditions = ",".join(dat[1].strip().split(",")[1:])
        self.typical_extreme_periods = ",".join(dat[2].strip().split(",")[1:])
        self.ground_temperatures = ",".join(dat[3].strip().split(",")[1:])
        self.holidays_daylight_savings = ",".join(dat[4].strip().split(",")[1:])
        self.comments_1 = ",".join(dat[5].strip().split(",")[1:])
        self.comments_2 = ",".join(dat[6].strip().split(",")[1:])

        # Read the data table
        df = pd.read_csv(StringIO("\n".join(dat[8:])), header=None)

        # Rename columns
        df.columns = [
            'year',
            'month',
            'day',
            'hour',
            'minute',
            'data_source_and_uncertainty_flags',
            'dry_bulb_temperature',
            'dew_point_temperature',
            'relative_humidity',
            'atmospheric_station_pressure',
            'extraterrestrial_horizontal_radiation',
            'extraterrestrial_direct_normal_radiation',
            'horizontal_infrared_radiation_intensity',
            'global_horizontal_radiation',
            'direct_normal_radiation',
            'diffuse_horizontal_radiation',
            'global_horizontal_illuminance',
            'direct_normal_illuminance',
            'diffuse_horizontal_illuminance',
            'zenith_luminance',
            'wind_direction',
            'wind_speed',
            'total_sky_cover',
            'opaque_sky_cover',
            'visibility',
            'ceiling_height',
            'present_weather_observation',
            'present_weather_codes',
            'precipitable_water',
            'aerosol_optical_depth',
            'snow_depth',
            'days_since_last_snowfall',
            'albedo',
            'liquid_precipitation_depth',
            'liquid_precipitation_quantity',
        ]

        # Create datetime index - using 2018 as base year (a Monday starting year without leap-day)
        self.index = self.index.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone)
        df.index = self.index

        # Drop date/time columns
        df.drop(columns=["year", "month", "day", "hour", "minute"], inplace=True)

        # Make loaded data accessible
        self.df = df
        self.data_source_and_uncertainty_flags = df.data_source_and_uncertainty_flags
        self.dry_bulb_temperature = df.dry_bulb_temperature
        self.dew_point_temperature = df.dew_point_temperature
        self.relative_humidity = df.relative_humidity
        self.atmospheric_station_pressure = df.atmospheric_station_pressure
        self.extraterrestrial_horizontal_radiation = df.extraterrestrial_horizontal_radiation
        self.extraterrestrial_direct_normal_radiation = df.extraterrestrial_direct_normal_radiation
        self.horizontal_infrared_radiation_intensity = df.horizontal_infrared_radiation_intensity
        self.global_horizontal_radiation = df.global_horizontal_radiation
        self.direct_normal_radiation = df.direct_normal_radiation
        self.diffuse_horizontal_radiation = df.diffuse_horizontal_radiation
        self.global_horizontal_illuminance = df.global_horizontal_illuminance
        self.direct_normal_illuminance = df.direct_normal_illuminance
        self.diffuse_horizontal_illuminance = df.diffuse_horizontal_illuminance
        self.zenith_luminance = df.zenith_luminance
        self.wind_direction = df.wind_direction
        self.wind_speed = df.wind_speed
        self.total_sky_cover = df.total_sky_cover
        self.opaque_sky_cover = df.opaque_sky_cover
        self.visibility = df.visibility
        self.ceiling_height = df.ceiling_height
        self.present_weather_observation = df.present_weather_observation
        self.present_weather_codes = df.present_weather_codes
        self.precipitable_water = df.precipitable_water
        self.aerosol_optical_depth = df.aerosol_optical_depth
        self.snow_depth = df.snow_depth
        self.days_since_last_snowfall = df.days_since_last_snowfall
        self.albedo = df.albedo
        self.liquid_precipitation_depth = df.liquid_precipitation_depth
        self.liquid_precipitation_quantity = df.liquid_precipitation_quantity

        # Calculate pedestrian height wind speed (at 1.5m above ground)
        self.pedestrian_wind_speed = pd.Series(name="pedestrian_wind_speed", index=self.index, data=[
            wind_speed_at_height(source_wind_speed=i, source_wind_height=10, target_wind_height=1.5) for i
            in self.wind_speed])
        self.df = pd.concat([self.df, self.pedestrian_wind_speed], axis=1)

        # Get the ground temperatures from the EPW file per month
        g_temps = {}
        for n, i in enumerate(list(chunk(self.ground_temperatures.split(",")[1:], 16))):
            g_temps[float(n)] = [float(j) for j in i[4:]]
        temp_ground_temperature = pd.DataFrame.from_dict(g_temps)
        temp_ground_temperature.index = pd.Series(index=self.index).resample("MS").mean().index
        temp_ground_temperature = pd.concat([pd.DataFrame(index=self.index), temp_ground_temperature], axis=1)
        temp_ground_temperature.columns = ["ground_temperature_0.5m", "ground_temperature_2m", "ground_temperature_4m"]
        temp_ground_temperature.iloc[-1, :] = temp_ground_temperature.iloc[0, :]  # Assign start temperature to last datetime
        temp_ground_temperature.interpolate(inplace=True)  # Fill in the gaps
        self.ground_temperature_1 = temp_ground_temperature["ground_temperature_0.5m"]
        self.ground_temperature_2 = temp_ground_temperature["ground_temperature_2m"]
        self.ground_temperature_3 = temp_ground_temperature["ground_temperature_4m"]
        self.df = pd.concat([self.df, temp_ground_temperature], axis=1)

        # Run the solar position calculations
        sol = sun_position(self.index, self.latitude, self.longitude)
        self.solar_apparent_zenith_angle = sol.solar_apparent_zenith_angle
        self.solar_zenith_angle = sol.solar_zenith_angle
        self.solar_apparent_elevation_angle = sol.solar_apparent_elevation_angle
        self.solar_elevation_angle = sol.solar_elevation_angle
        self.solar_azimuth_angle = sol.solar_azimuth_angle
        self.solar_equation_of_time = sol.solar_equation_of_time
        self.df = pd.concat([self.df, sol], axis=1)

        # Run the psychrometric calculations
        psych = psychrometric_calculations(self.dry_bulb_temperature, self.relative_humidity, self.atmospheric_station_pressure)
        psych.index = self.index
        self.humidity_ratio = psych.humidity_ratio
        self.wet_bulb_temperature = psych.wet_bulb_temperature
        self.partial_vapour_pressure_moist_air = psych.partial_vapour_pressure_moist_air
        self.enthalpy = psych.enthalpy
        self.specific_volume_moist_air = psych.specific_volume_moist_air
        self.degree_of_saturation = psych.degree_of_saturation
        self.df = pd.concat([self.df, psych], axis=1)

        if sky_matrix:

            # Specify sky patch sub-division method - Reinhart by default and hard-coded here
            self.reinhart = True

            # Create WEA file if it doesn't exist
            if self.wea_file is None:
                self.wea_file = self.to_wea()
                self.direct_sky_matrix_path = pathlib.Path(self.wea_file).with_suffix(".dirmtx")
                self.diffuse_sky_matrix_path = pathlib.Path(self.wea_file).with_suffix(".diffmtx")

            # Load pre-existing sky matrices if they exist. Currently no checks here for if the matrix matehces the weatherfile (other than
            if reuse_matrix:

                try:
                    self.direct_sky_matrix = load_sky_matrix(self.direct_sky_matrix_path, reinhart=self.reinhart)
                    self.diffuse_sky_matrix = load_sky_matrix(self.diffuse_sky_matrix_path, reinhart=self.reinhart)
                    self.total_sky_matrix = self.direct_sky_matrix + self.diffuse_sky_matrix
                    print("Direct and diffuse sky matrices loaded")
                except Exception as e:
                    raise ValueError("Looks like you haven't got a sky matrix to load - try creating one first!\n\n{0:}".format(e))
            else:
                # Generate sky matrices
                self.direct_sky_matrix, self.diffuse_sky_matrix, self.total_sky_matrix = sky_matrix_calculations(self.wea_file, reinhart=self.reinhart)

            # Set the other descriptors for the sky patches
            self.patch_centroids = REINHART_PATCH_CENTROIDS if self.reinhart else TREGENZA_PATCH_CENTROIDS
            self.patch_vectors = REINHART_PATCH_VECTORS if self.reinhart else TREGENZA_PATCH_VECTORS
            self.patch_count = REINHART_PATCH_COUNT if self.reinhart else TREGENZA_PATCH_COUNT

        # if openfield_mrt:
        #     if self.total_sky_matrix is None:
        #         self.sky_matrix(reinhart=True)
        #     self.mean_radiant_temperature_of()
        # if sa_mrt:
        #     self.mean_radiant_temperature_sa()

        return self

    def to_csv(self, file_path=None):
        """
        Write DataFrame containing hourly annual weather variables, and derived variables to a CSV file.

        Parameters
        ----------
        file_path : str
            The file path into which the annual hourly values will be written. If no value is passed, the output
            will be written to the same directory as the input weather-file.

        Returns
        -------
        str
            Path to the serialised CSV file.

        """
        if self.df is None:
            raise Exception('No data is available, try loading some first!')
        if file_path is None:
            file_path = pathlib.Path(self.file_path).with_suffix(".csv")
        self.df.to_csv(file_path)
        self.csv_file = str(file_path)
        print("CSV file created: {0:}".format(self.csv_file))
        return self.csv_file

    def to_wea(self, file_path=None):
        """
        Create WEA file for passing to Radiance programs.

        Parameters
        ----------
        file_path : str
            The file path into which the WEA formatted radiation values will be written. If no value is passed, the
            output will be written to the same directory as the input weather-file.

        Returns
        -------
        str
            Path to the serialised WEA file.

        """
        if (self.direct_normal_radiation is None) | (self.diffuse_horizontal_radiation is None):
            raise Exception('No radiation data is available, try loading some first!')
        if file_path is None:
            file_path = pathlib.Path(self.file_path).with_suffix(".wea")
        header = "place {0:}_{1:}\nlatitude {2:0.4f}\nlongitude {3:0.4f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.2f}\nweather_data_file_units 1".format(
            slugify(self.city), slugify(self.country), self.latitude, self.longitude, -self.time_zone / 15,
            self.elevation)
        values = []
        for n, dt in enumerate(self.index):
            values.append(
                "{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n],
                                                    self.diffuse_horizontal_radiation[n]))
        with open(file_path, "w") as f:
            f.write(header + "\n" + "\n".join(values) + "\n")
        self.wea_file = str(file_path)
        print("WEA file created: {0:}".format(self.wea_file))
        return self.wea_file

    def mean_radiant_temperature_openfield(self):
        mrt_openfield = mean_radiant_temperature_of(self)
        mrt_openfield.index = self.index
        mrt_openfield.name = "mean_radiant_temperature_openfield"
        self.mean_radiant_temperature_of = mrt_openfield
        self.df = pd.concat([self.df, self.mean_radiant_temperature_of], axis=1)

    def mean_radiant_temperature_solar_adjusted(self):
        mrt_solar_adjusted = mean_radiant_temperature_sa(self)
        mrt_solar_adjusted.index = self.index
        mrt_solar_adjusted.name = "mean_radiant_temperature_solar_adjusted"
        self.mean_radiant_temperature_sa = mrt_solar_adjusted
        self.df = pd.concat([self.df, self.mean_radiant_temperature_sa], axis=1)

    def utci_solar_adjusted(self):
        self.utci_sa = self.df.apply(
            lambda x: universal_thermal_climate_index(
                air_temperature=x.dry_bulb_temperature,
                mean_radiant_temperature=x.mean_radiant_temperature_solar_adjusted,
                air_velocity=x.pedestrian_wind_speed,
                relative_humidity=x.relative_humidity,
            ), axis=1
        )
        self.utci_sa.name = "utci_sa"
        self.df = pd.concat([self.df, self.utci_sa], axis=1)
        print("Universal thermal climate index calculations successful (SA)")

    def utci_openfield(self):
        self.utci_of = self.df.apply(
            lambda x: universal_thermal_climate_index(
                air_temperature=x.dry_bulb_temperature,
                mean_radiant_temperature=x.mean_radiant_temperature_openfield,
                air_velocity=x.pedestrian_wind_speed,
                relative_humidity=x.relative_humidity,
            ), axis=1
        )
        self.utci_of.name = "utci_of"
        self.df = pd.concat([self.df, self.utci_of], axis=1)
        print("Universal thermal climate index calculations successful (OF)")


# class Ground(object):
#
#     def __init__(self, material):
#         self.material = material
#         self.diffusivity = 0.001
#         self.thickness = 1
#         self.roughness = "Medium smooth"
#         self.emissivity = 0.8
#         self.absorptivity = 0.6
#         self.albedo = 0.2
#         self.surface_temperature = None
#         self.square_meter_radiation = None
#         self.k_value = None
#         self.convective_heat_transfer_coefficient = None
