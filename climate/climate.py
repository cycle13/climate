import platform
import subprocess

import pandas as pd
import numpy as np
import pathlib
from io import StringIO
from climate.helpers import slugify, chunk, angle_between, renamer
from pvlib.solarposition import get_solarposition
from psychrolib import SetUnitSystem, SI, CalcPsychrometricsFromRelHum

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from matplotlib import cm


class Weather(object):

    def __init__(self, file_path):

        # Metadata
        self.file_path = pathlib.Path(file_path).absolute()
        self.index = pd.date_range(
            start="2018-01-01 00:00:00",
            end="2019-01-01 00:00:00",
            freq="60T",
            closed="left"
        )

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

        # Path variables
        self.wea_file = None
        self.csv_file = None

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
        self.reinhart = None
        self.direct_sky_matrix = None
        self.diffuse_sky_matrix = None
        self.total_sky_matrix = None
        self.patch_centroids = None
        self.patch_vectors = None
        self.patch_count = None
        self.row_conversion_factor = None
        self.row_patches = None
        self.patch_conversion_factor = None
        self.radiation_rose_vectors = None
        self.radiation_rose_angles = None
        # Filtered rose values
        self.direct_sky_radiation_rose_values = None
        self.diffuse_sky_radiation_rose_values = None
        self.total_sky_radiation_rose_values = None

    # Methods below here for read/write
    def read(self, sun_position=False, psychrometrics=False, sky_matrix=False):
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
        self.index = self.index.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60))
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

        # Run the optional calculations
        if sun_position:
            self.sun_position()
        if psychrometrics:
            self.psychrometrics()
        if sky_matrix:
            self.sky_matrix()

        return self

    def to_csv(self, file_path=None):
        if self.df is None:
            raise Exception('No data is available, try loading some first!')
        if file_path is None:
            file_path = pathlib.Path(self.file_path).with_suffix(".csv")
        self.df.to_csv(file_path)
        self.csv_file = str(file_path)
        print("CSV file created: {0:}".format(self.csv_file))
        return self.csv_file

    def to_wea(self, file_path=None):
        if (self.direct_normal_radiation is None) | (self.diffuse_horizontal_radiation is None):
            raise Exception('No radiation data is available, try loading some first!')
        if file_path is None:
            file_path = pathlib.Path(self.file_path).with_suffix(".wea")
        header = "place {0:}_{1:}\nlatitude {2:0.4f}\nlongitude {3:0.4f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.2f}\nweather_data_file_units 1".format(slugify(self.city), slugify(self.country), self.latitude, self.longitude, -self.time_zone / 15, self.elevation)
        values = []
        for n, dt in enumerate(self.index):
            values.append("{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n], self.diffuse_horizontal_radiation[n]))
        with open(file_path, "w") as f:
            f.write(header + "\n" + "\n".join(values) + "\n")
        self.wea_file = str(file_path)
        print("WEA file created: {0:}".format(self.wea_file))
        return self.wea_file

    # Methods below here for derived psychrometrics/solar positioning

    def sun_position(self):
        solar_metrics = get_solarposition(self.index, self.latitude, self.longitude)
        solar_metrics.rename(columns={
            'apparent_zenith': 'solar_apparent_zenith_angle',
            'zenith': 'solar_zenith_angle',
            'apparent_elevation': 'solar_apparent_elevation_angle',
            'elevation': 'solar_elevation_angle',
            'azimuth': 'solar_azimuth_angle',
            'equation_of_time': 'solar_equation_of_time'
        }, inplace=True)
        self.solar_apparent_zenith_angle = solar_metrics.solar_apparent_zenith_angle
        self.solar_zenith_angle = solar_metrics.solar_zenith_angle
        self.solar_apparent_elevation_angle = solar_metrics.solar_apparent_elevation_angle
        self.solar_elevation_angle = solar_metrics.solar_elevation_angle
        self.solar_azimuth_angle = solar_metrics.solar_azimuth_angle
        self.solar_equation_of_time = solar_metrics.solar_equation_of_time
        self.df = pd.concat([self.df, solar_metrics], axis=1)
        print("Solar position calculations successful")
        return self

    def psychrometrics(self):
        SetUnitSystem(SI)
        psych_metrics = self.df.apply(lambda row: CalcPsychrometricsFromRelHum(row.dry_bulb_temperature, row.relative_humidity/100, row.atmospheric_station_pressure), axis=1).apply(pd.Series)
        psych_metrics.columns = ["humidity_ratio", "wet_bulb_temperature", "dew_point_temperature", "partial_vapour_pressure_moist_air", "enthalpy", "specific_volume_moist_air", "degree_of_saturation", ]
        self.dew_point_temperature = psych_metrics.dew_point_temperature
        self.humidity_ratio = psych_metrics.humidity_ratio
        self.wet_bulb_temperature = psych_metrics.wet_bulb_temperature
        self.partial_vapour_pressure_moist_air = psych_metrics.partial_vapour_pressure_moist_air
        self.enthalpy = psych_metrics.enthalpy
        self.specific_volume_moist_air = psych_metrics.specific_volume_moist_air
        self.degree_of_saturation = psych_metrics.degree_of_saturation
        self.df = pd.concat([self.df, psych_metrics], axis=1)
        print("Psychrometric calculations successful")
        return self

    def gendaymtx(self, direct=True):

        # Create output file path
        if direct:
            mtx_file = pathlib.Path(self.wea_file).with_suffix(".dirmtx")
        else:
            mtx_file = pathlib.Path(self.wea_file).with_suffix(".diffmtx")

        # Create run command
        if platform.system() != "Windows":
            cmd = '"/usr/local/radiance/bin/gendaymtx" -m {0:} -{1:} -O1 "{2:}" > "{3:}"'.format(2 if self.reinhart else 1,
                                                                                                 "d" if direct else "s",
                                                                                                 self.wea_file,
                                                                                                 mtx_file)
        else:
            cmd = '"C:/Radiance/bin/gendaymtx" -m {0:} -{1:} -O1 "{2:}" > "{3:}"'.format(2 if self.reinhart else 1,
                                                                                         "d" if direct else "s",
                                                                                         self.wea_file, mtx_file)

        # Run command
        subprocess.run(cmd, shell=True)

        # Load the resultant annual patch-value matrix
        sky_matrix = pd.read_csv(mtx_file, sep="\s+", skip_blank_lines=True, skiprows=8, header=None).values
        sky_matrix = np.sum(sky_matrix * np.array([0.265074126, 0.670114631, 0.064811243]), axis=1)
        sky_matrix = np.array(list(chunk(sky_matrix, 8760))[1:])
        sky_matrix = np.multiply(sky_matrix.T, self.patch_conversion_factor)

        return sky_matrix

    # Methods below here for sky-dome calculation

    def sky_matrix(self, reinhart=False):
        # Check if wea-file exists
        if self.wea_file is None:
            self.to_wea()

        self.reinhart = reinhart
        # Set sky-matrix patch metrics
        if reinhart:
            self.patch_centroids = np.array([
                [0, 0.997474, 0.054056],
                [0.104264, 0.99201, 0.054056],
                [0.207387, 0.975677, 0.054056],
                [0.308237, 0.948654, 0.054056],
                [0.405709, 0.911238, 0.054056],
                [0.498737, 0.863838, 0.054056],
                [0.586301, 0.806974, 0.054056],
                [0.667441, 0.741268, 0.054056],
                [0.741268, 0.667441, 0.054056],
                [0.806974, 0.586301, 0.054056],
                [0.863838, 0.498737, 0.054056],
                [0.911238, 0.405709, 0.054056],
                [0.948654, 0.308237, 0.054056],
                [0.975677, 0.207387, 0.054056],
                [0.99201, 0.104264, 0.054056],
                [0.997474, 0, 0.054056],
                [0.99201, -0.104264, 0.054056],
                [0.975677, -0.207387, 0.054056],
                [0.948654, -0.308237, 0.054056],
                [0.911238, -0.405709, 0.054056],
                [0.863838, -0.498737, 0.054056],
                [0.806974, -0.586301, 0.054056],
                [0.741268, -0.667441, 0.054056],
                [0.667441, -0.741268, 0.054056],
                [0.586301, -0.806974, 0.054056],
                [0.498737, -0.863838, 0.054056],
                [0.405709, -0.911238, 0.054056],
                [0.308237, -0.948654, 0.054056],
                [0.207387, -0.975677, 0.054056],
                [0.104264, -0.99201, 0.054056],
                [0, -0.997474, 0.054056],
                [-0.104264, -0.99201, 0.054056],
                [-0.207387, -0.975677, 0.054056],
                [-0.308237, -0.948654, 0.054056],
                [-0.405709, -0.911238, 0.054056],
                [-0.498737, -0.863838, 0.054056],
                [-0.586301, -0.806974, 0.054056],
                [-0.667441, -0.741268, 0.054056],
                [-0.741268, -0.667441, 0.054056],
                [-0.806974, -0.586301, 0.054056],
                [-0.863838, -0.498737, 0.054056],
                [-0.911238, -0.405709, 0.054056],
                [-0.948654, -0.308237, 0.054056],
                [-0.975677, -0.207387, 0.054056],
                [-0.99201, -0.104264, 0.054056],
                [-0.997474, 0, 0.054056],
                [-0.99201, 0.104264, 0.054056],
                [-0.975677, 0.207387, 0.054056],
                [-0.948654, 0.308237, 0.054056],
                [-0.911238, 0.405709, 0.054056],
                [-0.863838, 0.498737, 0.054056],
                [-0.806974, 0.586301, 0.054056],
                [-0.741268, 0.667441, 0.054056],
                [-0.667441, 0.741268, 0.054056],
                [-0.586301, 0.806974, 0.054056],
                [-0.498737, 0.863838, 0.054056],
                [-0.405709, 0.911238, 0.054056],
                [-0.308237, 0.948654, 0.054056],
                [-0.207387, 0.975677, 0.054056],
                [-0.104264, 0.99201, 0.054056],
                [0, 0.985803, 0.161535],
                [0.103044, 0.980403, 0.161535],
                [0.20496, 0.964261, 0.161535],
                [0.30463, 0.937554, 0.161535],
                [0.400962, 0.900576, 0.161535],
                [0.492901, 0.85373, 0.161535],
                [0.57944, 0.797531, 0.161535],
                [0.659631, 0.732594, 0.161535],
                [0.732594, 0.659631, 0.161535],
                [0.797531, 0.57944, 0.161535],
                [0.85373, 0.492901, 0.161535],
                [0.900576, 0.400962, 0.161535],
                [0.937554, 0.30463, 0.161535],
                [0.964261, 0.20496, 0.161535],
                [0.980403, 0.103044, 0.161535],
                [0.985803, 0, 0.161535],
                [0.980403, -0.103044, 0.161535],
                [0.964261, -0.20496, 0.161535],
                [0.937554, -0.30463, 0.161535],
                [0.900576, -0.400962, 0.161535],
                [0.85373, -0.492901, 0.161535],
                [0.797531, -0.57944, 0.161535],
                [0.732594, -0.659631, 0.161535],
                [0.659631, -0.732594, 0.161535],
                [0.57944, -0.797531, 0.161535],
                [0.492901, -0.85373, 0.161535],
                [0.400962, -0.900576, 0.161535],
                [0.30463, -0.937554, 0.161535],
                [0.20496, -0.964261, 0.161535],
                [0.103044, -0.980403, 0.161535],
                [0, -0.985803, 0.161535],
                [-0.103044, -0.980403, 0.161535],
                [-0.20496, -0.964261, 0.161535],
                [-0.30463, -0.937554, 0.161535],
                [-0.400962, -0.900576, 0.161535],
                [-0.492901, -0.85373, 0.161535],
                [-0.57944, -0.797531, 0.161535],
                [-0.659631, -0.732594, 0.161535],
                [-0.732594, -0.659631, 0.161535],
                [-0.797531, -0.57944, 0.161535],
                [-0.85373, -0.492901, 0.161535],
                [-0.900576, -0.400962, 0.161535],
                [-0.937554, -0.30463, 0.161535],
                [-0.964261, -0.20496, 0.161535],
                [-0.980403, -0.103044, 0.161535],
                [-0.985803, 0, 0.161535],
                [-0.980403, 0.103044, 0.161535],
                [-0.964261, 0.20496, 0.161535],
                [-0.937554, 0.30463, 0.161535],
                [-0.900576, 0.400962, 0.161535],
                [-0.85373, 0.492901, 0.161535],
                [-0.797531, 0.57944, 0.161535],
                [-0.732594, 0.659631, 0.161535],
                [-0.659631, 0.732594, 0.161535],
                [-0.57944, 0.797531, 0.161535],
                [-0.492901, 0.85373, 0.161535],
                [-0.400962, 0.900576, 0.161535],
                [-0.30463, 0.937554, 0.161535],
                [-0.20496, 0.964261, 0.161535],
                [-0.103044, 0.980403, 0.161535],
                [0, 0.962598, 0.26712],
                [0.100619, 0.957325, 0.26712],
                [0.200135, 0.941563, 0.26712],
                [0.297459, 0.915485, 0.26712],
                [0.391524, 0.879377, 0.26712],
                [0.481299, 0.833634, 0.26712],
                [0.565801, 0.778758, 0.26712],
                [0.644104, 0.71535, 0.26712],
                [0.71535, 0.644104, 0.26712],
                [0.778758, 0.565801, 0.26712],
                [0.833634, 0.481299, 0.26712],
                [0.879377, 0.391524, 0.26712],
                [0.915485, 0.297459, 0.26712],
                [0.941563, 0.200135, 0.26712],
                [0.957325, 0.100619, 0.26712],
                [0.962598, 0, 0.26712],
                [0.957325, -0.100619, 0.26712],
                [0.941563, -0.200135, 0.26712],
                [0.915485, -0.297459, 0.26712],
                [0.879377, -0.391524, 0.26712],
                [0.833634, -0.481299, 0.26712],
                [0.778758, -0.565801, 0.26712],
                [0.71535, -0.644104, 0.26712],
                [0.644104, -0.71535, 0.26712],
                [0.565801, -0.778758, 0.26712],
                [0.481299, -0.833634, 0.26712],
                [0.391524, -0.879377, 0.26712],
                [0.297459, -0.915485, 0.26712],
                [0.200135, -0.941563, 0.26712],
                [0.100619, -0.957325, 0.26712],
                [0, -0.962598, 0.26712],
                [-0.100619, -0.957325, 0.26712],
                [-0.200135, -0.941563, 0.26712],
                [-0.297459, -0.915485, 0.26712],
                [-0.391524, -0.879377, 0.26712],
                [-0.481299, -0.833634, 0.26712],
                [-0.565801, -0.778758, 0.26712],
                [-0.644104, -0.71535, 0.26712],
                [-0.71535, -0.644104, 0.26712],
                [-0.778758, -0.565801, 0.26712],
                [-0.833634, -0.481299, 0.26712],
                [-0.879377, -0.391524, 0.26712],
                [-0.915485, -0.297459, 0.26712],
                [-0.941563, -0.200135, 0.26712],
                [-0.957325, -0.100619, 0.26712],
                [-0.962598, 0, 0.26712],
                [-0.957325, 0.100619, 0.26712],
                [-0.941563, 0.200135, 0.26712],
                [-0.915485, 0.297459, 0.26712],
                [-0.879377, 0.391524, 0.26712],
                [-0.833634, 0.481299, 0.26712],
                [-0.778758, 0.565801, 0.26712],
                [-0.71535, 0.644104, 0.26712],
                [-0.644104, 0.71535, 0.26712],
                [-0.565801, 0.778758, 0.26712],
                [-0.481299, 0.833634, 0.26712],
                [-0.391524, 0.879377, 0.26712],
                [-0.297459, 0.915485, 0.26712],
                [-0.200135, 0.941563, 0.26712],
                [-0.100619, 0.957325, 0.26712],
                [0, 0.928133, 0.369573],
                [0.097016, 0.923049, 0.369573],
                [0.19297, 0.907851, 0.369573],
                [0.286809, 0.882707, 0.369573],
                [0.377506, 0.847892, 0.369573],
                [0.464066, 0.803787, 0.369573],
                [0.545543, 0.750875, 0.369573],
                [0.621042, 0.689737, 0.369573],
                [0.689737, 0.621042, 0.369573],
                [0.750875, 0.545543, 0.369573],
                [0.803787, 0.464066, 0.369573],
                [0.847892, 0.377506, 0.369573],
                [0.882707, 0.286809, 0.369573],
                [0.907851, 0.19297, 0.369573],
                [0.923049, 0.097016, 0.369573],
                [0.928133, 0, 0.369573],
                [0.923049, -0.097016, 0.369573],
                [0.907851, -0.19297, 0.369573],
                [0.882707, -0.286809, 0.369573],
                [0.847892, -0.377506, 0.369573],
                [0.803787, -0.464066, 0.369573],
                [0.750875, -0.545543, 0.369573],
                [0.689737, -0.621042, 0.369573],
                [0.621042, -0.689737, 0.369573],
                [0.545543, -0.750875, 0.369573],
                [0.464066, -0.803787, 0.369573],
                [0.377506, -0.847892, 0.369573],
                [0.286809, -0.882707, 0.369573],
                [0.19297, -0.907851, 0.369573],
                [0.097016, -0.923049, 0.369573],
                [0, -0.928133, 0.369573],
                [-0.097016, -0.923049, 0.369573],
                [-0.19297, -0.907851, 0.369573],
                [-0.286809, -0.882707, 0.369573],
                [-0.377506, -0.847892, 0.369573],
                [-0.464066, -0.803787, 0.369573],
                [-0.545543, -0.750875, 0.369573],
                [-0.621042, -0.689737, 0.369573],
                [-0.689737, -0.621042, 0.369573],
                [-0.750875, -0.545543, 0.369573],
                [-0.803787, -0.464066, 0.369573],
                [-0.847892, -0.377506, 0.369573],
                [-0.882707, -0.286809, 0.369573],
                [-0.907851, -0.19297, 0.369573],
                [-0.923049, -0.097016, 0.369573],
                [-0.928133, 0, 0.369573],
                [-0.923049, 0.097016, 0.369573],
                [-0.907851, 0.19297, 0.369573],
                [-0.882707, 0.286809, 0.369573],
                [-0.847892, 0.377506, 0.369573],
                [-0.803787, 0.464066, 0.369573],
                [-0.750875, 0.545543, 0.369573],
                [-0.689737, 0.621042, 0.369573],
                [-0.621042, 0.689737, 0.369573],
                [-0.545543, 0.750875, 0.369573],
                [-0.464066, 0.803787, 0.369573],
                [-0.377506, 0.847892, 0.369573],
                [-0.286809, 0.882707, 0.369573],
                [-0.19297, 0.907851, 0.369573],
                [-0.097016, 0.923049, 0.369573],
                [0, 0.882619, 0.467693],
                [0.115205, 0.875068, 0.467693],
                [0.228439, 0.852544, 0.467693],
                [0.337764, 0.815434, 0.467693],
                [0.441309, 0.76437, 0.467693],
                [0.537304, 0.700229, 0.467693],
                [0.624106, 0.624106, 0.467693],
                [0.700229, 0.537304, 0.467693],
                [0.76437, 0.441309, 0.467693],
                [0.815434, 0.337764, 0.467693],
                [0.852544, 0.228439, 0.467693],
                [0.875068, 0.115205, 0.467693],
                [0.882619, 0, 0.467693],
                [0.875068, -0.115205, 0.467693],
                [0.852544, -0.228439, 0.467693],
                [0.815434, -0.337764, 0.467693],
                [0.76437, -0.441309, 0.467693],
                [0.700229, -0.537304, 0.467693],
                [0.624106, -0.624106, 0.467693],
                [0.537304, -0.700229, 0.467693],
                [0.441309, -0.76437, 0.467693],
                [0.337764, -0.815434, 0.467693],
                [0.228439, -0.852544, 0.467693],
                [0.115205, -0.875068, 0.467693],
                [0, -0.882619, 0.467693],
                [-0.115205, -0.875068, 0.467693],
                [-0.228439, -0.852544, 0.467693],
                [-0.337764, -0.815434, 0.467693],
                [-0.441309, -0.76437, 0.467693],
                [-0.537304, -0.700229, 0.467693],
                [-0.624106, -0.624106, 0.467693],
                [-0.700229, -0.537304, 0.467693],
                [-0.76437, -0.441309, 0.467693],
                [-0.815434, -0.337764, 0.467693],
                [-0.852544, -0.228439, 0.467693],
                [-0.875068, -0.115205, 0.467693],
                [-0.882619, 0, 0.467693],
                [-0.875068, 0.115205, 0.467693],
                [-0.852544, 0.228439, 0.467693],
                [-0.815434, 0.337764, 0.467693],
                [-0.76437, 0.441309, 0.467693],
                [-0.700229, 0.537304, 0.467693],
                [-0.624106, 0.624106, 0.467693],
                [-0.537304, 0.700229, 0.467693],
                [-0.441309, 0.76437, 0.467693],
                [-0.337764, 0.815434, 0.467693],
                [-0.228439, 0.852544, 0.467693],
                [-0.115205, 0.875068, 0.467693],
                [0, 0.826997, 0.56033],
                [0.107945, 0.819922, 0.56033],
                [0.214043, 0.798818, 0.56033],
                [0.316478, 0.764045, 0.56033],
                [0.413498, 0.7162, 0.56033],
                [0.503444, 0.656101, 0.56033],
                [0.584775, 0.584775, 0.56033],
                [0.656101, 0.503444, 0.56033],
                [0.7162, 0.413498, 0.56033],
                [0.764045, 0.316478, 0.56033],
                [0.798818, 0.214043, 0.56033],
                [0.819922, 0.107945, 0.56033],
                [0.826997, 0, 0.56033],
                [0.819922, -0.107945, 0.56033],
                [0.798818, -0.214043, 0.56033],
                [0.764045, -0.316478, 0.56033],
                [0.7162, -0.413498, 0.56033],
                [0.656101, -0.503444, 0.56033],
                [0.584775, -0.584775, 0.56033],
                [0.503444, -0.656101, 0.56033],
                [0.413498, -0.7162, 0.56033],
                [0.316478, -0.764045, 0.56033],
                [0.214043, -0.798818, 0.56033],
                [0.107945, -0.819922, 0.56033],
                [0, -0.826997, 0.56033],
                [-0.107945, -0.819922, 0.56033],
                [-0.214043, -0.798818, 0.56033],
                [-0.316478, -0.764045, 0.56033],
                [-0.413498, -0.7162, 0.56033],
                [-0.503444, -0.656101, 0.56033],
                [-0.584775, -0.584775, 0.56033],
                [-0.656101, -0.503444, 0.56033],
                [-0.7162, -0.413498, 0.56033],
                [-0.764045, -0.316478, 0.56033],
                [-0.798818, -0.214043, 0.56033],
                [-0.819922, -0.107945, 0.56033],
                [-0.826997, 0, 0.56033],
                [-0.819922, 0.107945, 0.56033],
                [-0.798818, 0.214043, 0.56033],
                [-0.764045, 0.316478, 0.56033],
                [-0.7162, 0.413498, 0.56033],
                [-0.656101, 0.503444, 0.56033],
                [-0.584775, 0.584775, 0.56033],
                [-0.503444, 0.656101, 0.56033],
                [-0.413498, 0.7162, 0.56033],
                [-0.316478, 0.764045, 0.56033],
                [-0.214043, 0.798818, 0.56033],
                [-0.107945, 0.819922, 0.56033],
                [0, 0.76172, 0.646397],
                [0.099424, 0.755203, 0.646397],
                [0.197148, 0.735765, 0.646397],
                [0.291497, 0.703737, 0.646397],
                [0.38086, 0.659669, 0.646397],
                [0.463706, 0.604313, 0.646397],
                [0.538617, 0.538617, 0.646397],
                [0.604313, 0.463706, 0.646397],
                [0.659669, 0.38086, 0.646397],
                [0.703737, 0.291497, 0.646397],
                [0.735765, 0.197148, 0.646397],
                [0.755203, 0.099424, 0.646397],
                [0.76172, 0, 0.646397],
                [0.755203, -0.099424, 0.646397],
                [0.735765, -0.197148, 0.646397],
                [0.703737, -0.291497, 0.646397],
                [0.659669, -0.38086, 0.646397],
                [0.604313, -0.463706, 0.646397],
                [0.538617, -0.538617, 0.646397],
                [0.463706, -0.604313, 0.646397],
                [0.38086, -0.659669, 0.646397],
                [0.291497, -0.703737, 0.646397],
                [0.197148, -0.735765, 0.646397],
                [0.099424, -0.755203, 0.646397],
                [0, -0.76172, 0.646397],
                [-0.099424, -0.755203, 0.646397],
                [-0.197148, -0.735765, 0.646397],
                [-0.291497, -0.703737, 0.646397],
                [-0.38086, -0.659669, 0.646397],
                [-0.463706, -0.604313, 0.646397],
                [-0.538617, -0.538617, 0.646397],
                [-0.604313, -0.463706, 0.646397],
                [-0.659669, -0.38086, 0.646397],
                [-0.703737, -0.291497, 0.646397],
                [-0.735765, -0.197148, 0.646397],
                [-0.755203, -0.099424, 0.646397],
                [-0.76172, 0, 0.646397],
                [-0.755203, 0.099424, 0.646397],
                [-0.735765, 0.197148, 0.646397],
                [-0.703737, 0.291497, 0.646397],
                [-0.659669, 0.38086, 0.646397],
                [-0.604313, 0.463706, 0.646397],
                [-0.538617, 0.538617, 0.646397],
                [-0.463706, 0.604313, 0.646397],
                [-0.38086, 0.659669, 0.646397],
                [-0.291497, 0.703737, 0.646397],
                [-0.197148, 0.735765, 0.646397],
                [-0.099424, 0.755203, 0.646397],
                [0, 0.687518, 0.724886],
                [0.089739, 0.681637, 0.724886],
                [0.177943, 0.664092, 0.724886],
                [0.263102, 0.635184, 0.724886],
                [0.343759, 0.595408, 0.724886],
                [0.418535, 0.545445, 0.724886],
                [0.486149, 0.486149, 0.724886],
                [0.545445, 0.418535, 0.724886],
                [0.595408, 0.343759, 0.724886],
                [0.635184, 0.263102, 0.724886],
                [0.664092, 0.177943, 0.724886],
                [0.681637, 0.089739, 0.724886],
                [0.687518, 0, 0.724886],
                [0.681637, -0.089739, 0.724886],
                [0.664092, -0.177943, 0.724886],
                [0.635184, -0.263102, 0.724886],
                [0.595408, -0.343759, 0.724886],
                [0.545445, -0.418535, 0.724886],
                [0.486149, -0.486149, 0.724886],
                [0.418535, -0.545445, 0.724886],
                [0.343759, -0.595408, 0.724886],
                [0.263102, -0.635184, 0.724886],
                [0.177943, -0.664092, 0.724886],
                [0.089739, -0.681637, 0.724886],
                [0, -0.687518, 0.724886],
                [-0.089739, -0.681637, 0.724886],
                [-0.177943, -0.664092, 0.724886],
                [-0.263102, -0.635184, 0.724886],
                [-0.343759, -0.595408, 0.724886],
                [-0.418535, -0.545445, 0.724886],
                [-0.486149, -0.486149, 0.724886],
                [-0.545445, -0.418535, 0.724886],
                [-0.595408, -0.343759, 0.724886],
                [-0.635184, -0.263102, 0.724886],
                [-0.664092, -0.177943, 0.724886],
                [-0.681637, -0.089739, 0.724886],
                [-0.687518, 0, 0.724886],
                [-0.681637, 0.089739, 0.724886],
                [-0.664092, 0.177943, 0.724886],
                [-0.635184, 0.263102, 0.724886],
                [-0.595408, 0.343759, 0.724886],
                [-0.545445, 0.418535, 0.724886],
                [-0.486149, 0.486149, 0.724886],
                [-0.418535, 0.545445, 0.724886],
                [-0.343759, 0.595408, 0.724886],
                [-0.263102, 0.635184, 0.724886],
                [-0.177943, 0.664092, 0.724886],
                [-0.089739, 0.681637, 0.724886],
                [0, 0.605073, 0.794877],
                [0.10507, 0.59588, 0.794877],
                [0.206947, 0.568582, 0.794877],
                [0.302536, 0.524008, 0.794877],
                [0.388933, 0.463513, 0.794877],
                [0.463513, 0.388933, 0.794877],
                [0.524008, 0.302536, 0.794877],
                [0.568582, 0.206947, 0.794877],
                [0.59588, 0.10507, 0.794877],
                [0.605073, 0, 0.794877],
                [0.59588, -0.10507, 0.794877],
                [0.568582, -0.206947, 0.794877],
                [0.524008, -0.302536, 0.794877],
                [0.463513, -0.388933, 0.794877],
                [0.388933, -0.463513, 0.794877],
                [0.302536, -0.524008, 0.794877],
                [0.206947, -0.568582, 0.794877],
                [0.10507, -0.59588, 0.794877],
                [0, -0.605073, 0.794877],
                [-0.10507, -0.59588, 0.794877],
                [-0.206947, -0.568582, 0.794877],
                [-0.302536, -0.524008, 0.794877],
                [-0.388933, -0.463513, 0.794877],
                [-0.463513, -0.388933, 0.794877],
                [-0.524008, -0.302536, 0.794877],
                [-0.568582, -0.206947, 0.794877],
                [-0.59588, -0.10507, 0.794877],
                [-0.605073, 0, 0.794877],
                [-0.59588, 0.10507, 0.794877],
                [-0.568582, 0.206947, 0.794877],
                [-0.524008, 0.302536, 0.794877],
                [-0.463513, 0.388933, 0.794877],
                [-0.388933, 0.463513, 0.794877],
                [-0.302536, 0.524008, 0.794877],
                [-0.206947, 0.568582, 0.794877],
                [-0.10507, 0.59588, 0.794877],
                [0, 0.515987, 0.855548],
                [0.0896, 0.508148, 0.855548],
                [0.176478, 0.484869, 0.855548],
                [0.257993, 0.446858, 0.855548],
                [0.33167, 0.395269, 0.855548],
                [0.395269, 0.33167, 0.855548],
                [0.446858, 0.257993, 0.855548],
                [0.484869, 0.176478, 0.855548],
                [0.508148, 0.0896, 0.855548],
                [0.515987, 0, 0.855548],
                [0.508148, -0.0896, 0.855548],
                [0.484869, -0.176478, 0.855548],
                [0.446858, -0.257993, 0.855548],
                [0.395269, -0.33167, 0.855548],
                [0.33167, -0.395269, 0.855548],
                [0.257993, -0.446858, 0.855548],
                [0.176478, -0.484869, 0.855548],
                [0.0896, -0.508148, 0.855548],
                [0, -0.515987, 0.855548],
                [-0.0896, -0.508148, 0.855548],
                [-0.176478, -0.484869, 0.855548],
                [-0.257993, -0.446858, 0.855548],
                [-0.33167, -0.395269, 0.855548],
                [-0.395269, -0.33167, 0.855548],
                [-0.446858, -0.257993, 0.855548],
                [-0.484869, -0.176478, 0.855548],
                [-0.508148, -0.0896, 0.855548],
                [-0.515987, 0, 0.855548],
                [-0.508148, 0.0896, 0.855548],
                [-0.484869, 0.176478, 0.855548],
                [-0.446858, 0.257993, 0.855548],
                [-0.395269, 0.33167, 0.855548],
                [-0.33167, 0.395269, 0.855548],
                [-0.257993, 0.446858, 0.855548],
                [-0.176478, 0.484869, 0.855548],
                [-0.0896, 0.508148, 0.855548],
                [0, 0.420336, 0.906189],
                [0.108791, 0.406013, 0.906189],
                [0.210168, 0.364022, 0.906189],
                [0.297222, 0.297222, 0.906189],
                [0.364022, 0.210168, 0.906189],
                [0.406013, 0.108791, 0.906189],
                [0.420336, 0, 0.906189],
                [0.406013, -0.108791, 0.906189],
                [0.364022, -0.210168, 0.906189],
                [0.297222, -0.297222, 0.906189],
                [0.210168, -0.364022, 0.906189],
                [0.108791, -0.406013, 0.906189],
                [0, -0.420336, 0.906189],
                [-0.108791, -0.406013, 0.906189],
                [-0.210168, -0.364022, 0.906189],
                [-0.297222, -0.297222, 0.906189],
                [-0.364022, -0.210168, 0.906189],
                [-0.406013, -0.108791, 0.906189],
                [-0.420336, 0, 0.906189],
                [-0.406013, 0.108791, 0.906189],
                [-0.364022, 0.210168, 0.906189],
                [-0.297222, 0.297222, 0.906189],
                [-0.210168, 0.364022, 0.906189],
                [-0.108791, 0.406013, 0.906189],
                [0, 0.320929, 0.946205],
                [0.083063, 0.309994, 0.946205],
                [0.160465, 0.277933, 0.946205],
                [0.226931, 0.226931, 0.946205],
                [0.277933, 0.160465, 0.946205],
                [0.309994, 0.083063, 0.946205],
                [0.320929, 0, 0.946205],
                [0.309994, -0.083063, 0.946205],
                [0.277933, -0.160465, 0.946205],
                [0.226931, -0.226931, 0.946205],
                [0.160465, -0.277933, 0.946205],
                [0.083063, -0.309994, 0.946205],
                [0, -0.320929, 0.946205],
                [-0.083063, -0.309994, 0.946205],
                [-0.160465, -0.277933, 0.946205],
                [-0.226931, -0.226931, 0.946205],
                [-0.277933, -0.160465, 0.946205],
                [-0.309994, -0.083063, 0.946205],
                [-0.320929, 0, 0.946205],
                [-0.309994, 0.083063, 0.946205],
                [-0.277933, 0.160465, 0.946205],
                [-0.226931, 0.226931, 0.946205],
                [-0.160465, 0.277933, 0.946205],
                [-0.083063, 0.309994, 0.946205],
                [0, 0.216676, 0.975129],
                [0.108338, 0.187647, 0.975129],
                [0.187647, 0.108338, 0.975129],
                [0.216676, 0, 0.975129],
                [0.187647, -0.108338, 0.975129],
                [0.108338, -0.187647, 0.975129],
                [0, -0.216676, 0.975129],
                [-0.108338, -0.187647, 0.975129],
                [-0.187647, -0.108338, 0.975129],
                [-0.216676, 0, 0.975129],
                [-0.187647, 0.108338, 0.975129],
                [-0.108338, 0.187647, 0.975129],
                [0, 0.115625, 0.992619],
                [0.057812, 0.100134, 0.992619],
                [0.100134, 0.057812, 0.992619],
                [0.115625, 0, 0.992619],
                [0.100134, -0.057812, 0.992619],
                [0.057812, -0.100134, 0.992619],
                [0, -0.115625, 0.992619],
                [-0.057812, -0.100134, 0.992619],
                [-0.100134, -0.057812, 0.992619],
                [-0.115625, 0, 0.992619],
                [-0.100134, 0.057812, 0.992619],
                [-0.057812, 0.100134, 0.992619],
                [0, 0, 0.999251],
            ])
            self.patch_vectors = np.array([
                [-0.0000000, 0.9985349, 0.0541114],
                [0.1043753, 0.9930648, 0.0541114],
                [0.2076071, 0.9767145, 0.0541114],
                [0.3085643, 0.9496631, 0.0541114],
                [0.4061407, 0.9122070, 0.0541114],
                [0.4992675, 0.8647566, 0.0541114],
                [0.5869241, 0.8078317, 0.0541114],
                [0.6681503, 0.7420561, 0.0541114],
                [0.7420561, 0.6681503, 0.0541114],
                [0.8078317, 0.5869241, 0.0541114],
                [0.8647566, 0.4992675, 0.0541114],
                [0.9122070, 0.4061407, 0.0541114],
                [0.9496631, 0.3085643, 0.0541114],
                [0.9767145, 0.2076071, 0.0541114],
                [0.9930648, 0.1043753, 0.0541114],
                [0.9985349, -0.0000000, 0.0541114],
                [0.9930648, -0.1043753, 0.0541114],
                [0.9767145, -0.2076071, 0.0541114],
                [0.9496631, -0.3085643, 0.0541114],
                [0.9122070, -0.4061407, 0.0541114],
                [0.8647566, -0.4992675, 0.0541114],
                [0.8078317, -0.5869241, 0.0541114],
                [0.7420561, -0.6681503, 0.0541114],
                [0.6681503, -0.7420561, 0.0541114],
                [0.5869241, -0.8078317, 0.0541114],
                [0.4992675, -0.8647566, 0.0541114],
                [0.4061407, -0.9122070, 0.0541114],
                [0.3085643, -0.9496631, 0.0541114],
                [0.2076071, -0.9767145, 0.0541114],
                [0.1043753, -0.9930648, 0.0541114],
                [-0.0000000, -0.9985349, 0.0541114],
                [-0.1043753, -0.9930648, 0.0541114],
                [-0.2076071, -0.9767145, 0.0541114],
                [-0.3085643, -0.9496631, 0.0541114],
                [-0.4061407, -0.9122070, 0.0541114],
                [-0.4992675, -0.8647566, 0.0541114],
                [-0.5869241, -0.8078317, 0.0541114],
                [-0.6681503, -0.7420561, 0.0541114],
                [-0.7420561, -0.6681503, 0.0541114],
                [-0.8078317, -0.5869241, 0.0541114],
                [-0.8647566, -0.4992675, 0.0541114],
                [-0.9122070, -0.4061407, 0.0541114],
                [-0.9496631, -0.3085643, 0.0541114],
                [-0.9767145, -0.2076071, 0.0541114],
                [-0.9930648, -0.1043753, 0.0541114],
                [-0.9985349, 0.0000000, 0.0541114],
                [-0.9930648, 0.1043753, 0.0541114],
                [-0.9767145, 0.2076071, 0.0541114],
                [-0.9496631, 0.3085643, 0.0541114],
                [-0.9122070, 0.4061407, 0.0541114],
                [-0.8647566, 0.4992675, 0.0541114],
                [-0.8078317, 0.5869241, 0.0541114],
                [-0.7420561, 0.6681503, 0.0541114],
                [-0.6681503, 0.7420561, 0.0541114],
                [-0.5869241, 0.8078317, 0.0541114],
                [-0.4992675, 0.8647566, 0.0541114],
                [-0.4061407, 0.9122070, 0.0541114],
                [-0.3085643, 0.9496631, 0.0541114],
                [-0.2076071, 0.9767145, 0.0541114],
                [-0.1043753, 0.9930648, 0.0541114],
                [0.0000000, 0.9868403, 0.1616979],
                [0.1031529, 0.9814343, 0.1616979],
                [0.2051756, 0.9652755, 0.1616979],
                [0.3049504, 0.9385409, 0.1616979],
                [0.4013841, 0.9015235, 0.1616979],
                [0.4934202, 0.8546288, 0.1616979],
                [0.5800502, 0.7983706, 0.1616979],
                [0.6603250, 0.7333653, 0.1616979],
                [0.7333653, 0.6603250, 0.1616979],
                [0.7983706, 0.5800502, 0.1616979],
                [0.8546288, 0.4934202, 0.1616979],
                [0.9015235, 0.4013841, 0.1616979],
                [0.9385409, 0.3049504, 0.1616979],
                [0.9652755, 0.2051756, 0.1616979],
                [0.9814343, 0.1031529, 0.1616979],
                [0.9868403, 0.0000000, 0.1616979],
                [0.9814343, -0.1031529, 0.1616979],
                [0.9652755, -0.2051756, 0.1616979],
                [0.9385409, -0.3049504, 0.1616979],
                [0.9015235, -0.4013841, 0.1616979],
                [0.8546288, -0.4934202, 0.1616979],
                [0.7983706, -0.5800502, 0.1616979],
                [0.7333653, -0.6603250, 0.1616979],
                [0.6603250, -0.7333653, 0.1616979],
                [0.5800502, -0.7983706, 0.1616979],
                [0.4934202, -0.8546288, 0.1616979],
                [0.4013841, -0.9015235, 0.1616979],
                [0.3049504, -0.9385409, 0.1616979],
                [0.2051756, -0.9652755, 0.1616979],
                [0.1031529, -0.9814343, 0.1616979],
                [0.0000000, -0.9868403, 0.1616979],
                [-0.1031529, -0.9814343, 0.1616979],
                [-0.2051756, -0.9652755, 0.1616979],
                [-0.3049504, -0.9385409, 0.1616979],
                [-0.4013841, -0.9015235, 0.1616979],
                [-0.4934202, -0.8546288, 0.1616979],
                [-0.5800502, -0.7983706, 0.1616979],
                [-0.6603250, -0.7333653, 0.1616979],
                [-0.7333653, -0.6603250, 0.1616979],
                [-0.7983706, -0.5800502, 0.1616979],
                [-0.8546288, -0.4934202, 0.1616979],
                [-0.9015235, -0.4013841, 0.1616979],
                [-0.9385409, -0.3049504, 0.1616979],
                [-0.9652755, -0.2051756, 0.1616979],
                [-0.9814343, -0.1031529, 0.1616979],
                [-0.9868403, 0.0000000, 0.1616979],
                [-0.9814343, 0.1031529, 0.1616979],
                [-0.9652755, 0.2051756, 0.1616979],
                [-0.9385409, 0.3049504, 0.1616979],
                [-0.9015235, 0.4013841, 0.1616979],
                [-0.8546288, 0.4934202, 0.1616979],
                [-0.7983706, 0.5800502, 0.1616979],
                [-0.7333653, 0.6603250, 0.1616979],
                [-0.6603250, 0.7333653, 0.1616979],
                [-0.5800502, 0.7983706, 0.1616979],
                [-0.4934202, 0.8546288, 0.1616979],
                [-0.4013841, 0.9015235, 0.1616979],
                [-0.3049504, 0.9385409, 0.1616979],
                [-0.2051756, 0.9652755, 0.1616979],
                [-0.1031529, 0.9814343, 0.1616979],
                [0.0000000, 0.9635902, 0.2673836],
                [0.1007226, 0.9583115, 0.2673836],
                [0.2003417, 0.9425334, 0.2673836],
                [0.2977657, 0.9164287, 0.2673836],
                [0.3919274, 0.8802834, 0.2673836],
                [0.4817951, 0.8344936, 0.2673836],
                [0.5663841, 0.7795608, 0.2673836],
                [0.6447677, 0.7160871, 0.2673836],
                [0.7160871, 0.6447677, 0.2673836],
                [0.7795608, 0.5663841, 0.2673836],
                [0.8344936, 0.4817951, 0.2673836],
                [0.8802834, 0.3919274, 0.2673836],
                [0.9164287, 0.2977657, 0.2673836],
                [0.9425334, 0.2003417, 0.2673836],
                [0.9583115, 0.1007226, 0.2673836],
                [0.9635902, 0.0000000, 0.2673836],
                [0.9583115, -0.1007226, 0.2673836],
                [0.9425334, -0.2003417, 0.2673836],
                [0.9164287, -0.2977657, 0.2673836],
                [0.8802834, -0.3919274, 0.2673836],
                [0.8344936, -0.4817951, 0.2673836],
                [0.7795608, -0.5663841, 0.2673836],
                [0.7160871, -0.6447677, 0.2673836],
                [0.6447677, -0.7160871, 0.2673836],
                [0.5663841, -0.7795608, 0.2673836],
                [0.4817951, -0.8344936, 0.2673836],
                [0.3919274, -0.8802834, 0.2673836],
                [0.2977657, -0.9164287, 0.2673836],
                [0.2003417, -0.9425334, 0.2673836],
                [0.1007226, -0.9583115, 0.2673836],
                [0.0000000, -0.9635902, 0.2673836],
                [-0.1007226, -0.9583115, 0.2673836],
                [-0.2003417, -0.9425334, 0.2673836],
                [-0.2977657, -0.9164287, 0.2673836],
                [-0.3919274, -0.8802834, 0.2673836],
                [-0.4817951, -0.8344936, 0.2673836],
                [-0.5663841, -0.7795608, 0.2673836],
                [-0.6447677, -0.7160871, 0.2673836],
                [-0.7160871, -0.6447677, 0.2673836],
                [-0.7795608, -0.5663841, 0.2673836],
                [-0.8344936, -0.4817951, 0.2673836],
                [-0.8802834, -0.3919274, 0.2673836],
                [-0.9164287, -0.2977657, 0.2673836],
                [-0.9425334, -0.2003417, 0.2673836],
                [-0.9583115, -0.1007226, 0.2673836],
                [-0.9635902, 0.0000000, 0.2673836],
                [-0.9583115, 0.1007226, 0.2673836],
                [-0.9425334, 0.2003417, 0.2673836],
                [-0.9164287, 0.2977657, 0.2673836],
                [-0.8802834, 0.3919274, 0.2673836],
                [-0.8344936, 0.4817951, 0.2673836],
                [-0.7795608, 0.5663841, 0.2673836],
                [-0.7160871, 0.6447677, 0.2673836],
                [-0.6447677, 0.7160871, 0.2673836],
                [-0.5663841, 0.7795608, 0.2673836],
                [-0.4817951, 0.8344936, 0.2673836],
                [-0.3919274, 0.8802834, 0.2673836],
                [-0.2977657, 0.9164287, 0.2673836],
                [-0.2003417, 0.9425334, 0.2673836],
                [-0.1007226, 0.9583115, 0.2673836],
                [0.0000000, 0.9290610, 0.3699264],
                [0.0971133, 0.9239716, 0.3699264],
                [0.1931627, 0.9087588, 0.3699264],
                [0.2870957, 0.8835896, 0.3699264],
                [0.3778832, 0.8487395, 0.3699264],
                [0.4645305, 0.8045905, 0.3699264],
                [0.5460884, 0.7516262, 0.3699264],
                [0.6216632, 0.6904269, 0.3699264],
                [0.6904269, 0.6216632, 0.3699264],
                [0.7516262, 0.5460884, 0.3699264],
                [0.8045905, 0.4645305, 0.3699264],
                [0.8487395, 0.3778832, 0.3699264],
                [0.8835896, 0.2870957, 0.3699264],
                [0.9087588, 0.1931627, 0.3699264],
                [0.9239716, 0.0971133, 0.3699264],
                [0.9290610, 0.0000000, 0.3699264],
                [0.9239716, -0.0971133, 0.3699264],
                [0.9087588, -0.1931627, 0.3699264],
                [0.8835896, -0.2870957, 0.3699264],
                [0.8487395, -0.3778832, 0.3699264],
                [0.8045905, -0.4645305, 0.3699264],
                [0.7516262, -0.5460884, 0.3699264],
                [0.6904269, -0.6216632, 0.3699264],
                [0.6216632, -0.6904269, 0.3699264],
                [0.5460884, -0.7516262, 0.3699264],
                [0.4645305, -0.8045905, 0.3699264],
                [0.3778832, -0.8487395, 0.3699264],
                [0.2870957, -0.8835896, 0.3699264],
                [0.1931627, -0.9087588, 0.3699264],
                [0.0971133, -0.9239716, 0.3699264],
                [0.0000000, -0.9290610, 0.3699264],
                [-0.0971133, -0.9239716, 0.3699264],
                [-0.1931627, -0.9087588, 0.3699264],
                [-0.2870957, -0.8835896, 0.3699264],
                [-0.3778832, -0.8487395, 0.3699264],
                [-0.4645305, -0.8045905, 0.3699264],
                [-0.5460884, -0.7516262, 0.3699264],
                [-0.6216632, -0.6904269, 0.3699264],
                [-0.6904269, -0.6216632, 0.3699264],
                [-0.7516262, -0.5460884, 0.3699264],
                [-0.8045905, -0.4645305, 0.3699264],
                [-0.8487395, -0.3778832, 0.3699264],
                [-0.8835896, -0.2870957, 0.3699264],
                [-0.9087588, -0.1931627, 0.3699264],
                [-0.9239716, -0.0971133, 0.3699264],
                [-0.9290610, 0.0000000, 0.3699264],
                [-0.9239716, 0.0971133, 0.3699264],
                [-0.9087588, 0.1931627, 0.3699264],
                [-0.8835896, 0.2870957, 0.3699264],
                [-0.8487395, 0.3778832, 0.3699264],
                [-0.8045905, 0.4645305, 0.3699264],
                [-0.7516262, 0.5460884, 0.3699264],
                [-0.6904269, 0.6216632, 0.3699264],
                [-0.6216632, 0.6904269, 0.3699264],
                [-0.5460884, 0.7516262, 0.3699264],
                [-0.4645305, 0.8045905, 0.3699264],
                [-0.3778832, 0.8487395, 0.3699264],
                [-0.2870957, 0.8835896, 0.3699264],
                [-0.1931627, 0.9087588, 0.3699264],
                [-0.0971133, 0.9239716, 0.3699264],
                [0.0000000, 0.8836123, 0.4682192],
                [0.1153345, 0.8760529, 0.4682192],
                [0.2286957, 0.8535040, 0.4682193],
                [0.3381438, 0.8163513, 0.4682192],
                [0.4418062, 0.7652307, 0.4682192],
                [0.5379091, 0.7010168, 0.4682192],
                [0.6248083, 0.6248083, 0.4682193],
                [0.7010168, 0.5379091, 0.4682192],
                [0.7652307, 0.4418062, 0.4682192],
                [0.8163513, 0.3381438, 0.4682192],
                [0.8535040, 0.2286957, 0.4682193],
                [0.8760529, 0.1153345, 0.4682192],
                [0.8836123, 0.0000000, 0.4682192],
                [0.8760529, -0.1153345, 0.4682192],
                [0.8535040, -0.2286957, 0.4682193],
                [0.8163513, -0.3381438, 0.4682192],
                [0.7652307, -0.4418062, 0.4682192],
                [0.7010168, -0.5379091, 0.4682192],
                [0.6248083, -0.6248083, 0.4682193],
                [0.5379091, -0.7010168, 0.4682192],
                [0.4418062, -0.7652307, 0.4682192],
                [0.3381438, -0.8163513, 0.4682192],
                [0.2286957, -0.8535040, 0.4682193],
                [0.1153345, -0.8760529, 0.4682192],
                [0.0000000, -0.8836123, 0.4682192],
                [-0.1153345, -0.8760529, 0.4682192],
                [-0.2286957, -0.8535040, 0.4682193],
                [-0.3381438, -0.8163513, 0.4682192],
                [-0.4418062, -0.7652307, 0.4682192],
                [-0.5379091, -0.7010168, 0.4682192],
                [-0.6248083, -0.6248083, 0.4682193],
                [-0.7010168, -0.5379091, 0.4682192],
                [-0.7652307, -0.4418062, 0.4682192],
                [-0.8163513, -0.3381438, 0.4682192],
                [-0.8535040, -0.2286957, 0.4682193],
                [-0.8760529, -0.1153345, 0.4682192],
                [-0.8836123, 0.0000000, 0.4682192],
                [-0.8760529, 0.1153345, 0.4682192],
                [-0.8535040, 0.2286957, 0.4682193],
                [-0.8163513, 0.3381438, 0.4682192],
                [-0.7652307, 0.4418062, 0.4682192],
                [-0.7010168, 0.5379091, 0.4682192],
                [-0.6248083, 0.6248083, 0.4682193],
                [-0.5379091, 0.7010168, 0.4682192],
                [-0.4418062, 0.7652307, 0.4682192],
                [-0.3381438, 0.8163513, 0.4682192],
                [-0.2286957, 0.8535040, 0.4682193],
                [-0.1153345, 0.8760529, 0.4682192],
                [0.0000000, 0.8278694, 0.5609210],
                [0.1080586, 0.8207868, 0.5609210],
                [0.2142684, 0.7996604, 0.5609210],
                [0.3168119, 0.7648516, 0.5609210],
                [0.4139347, 0.7169559, 0.5609210],
                [0.5039749, 0.6567929, 0.5609210],
                [0.5853920, 0.5853920, 0.5609210],
                [0.6567929, 0.5039749, 0.5609210],
                [0.7169559, 0.4139347, 0.5609210],
                [0.7648516, 0.3168119, 0.5609210],
                [0.7996604, 0.2142684, 0.5609210],
                [0.8207868, 0.1080586, 0.5609210],
                [0.8278694, 0.0000000, 0.5609210],
                [0.8207868, -0.1080586, 0.5609210],
                [0.7996604, -0.2142684, 0.5609210],
                [0.7648516, -0.3168119, 0.5609210],
                [0.7169559, -0.4139347, 0.5609210],
                [0.6567929, -0.5039749, 0.5609210],
                [0.5853920, -0.5853920, 0.5609210],
                [0.5039749, -0.6567929, 0.5609210],
                [0.4139347, -0.7169559, 0.5609210],
                [0.3168119, -0.7648516, 0.5609210],
                [0.2142684, -0.7996604, 0.5609210],
                [0.1080586, -0.8207868, 0.5609210],
                [0.0000000, -0.8278694, 0.5609210],
                [-0.1080586, -0.8207868, 0.5609210],
                [-0.2142684, -0.7996604, 0.5609210],
                [-0.3168119, -0.7648516, 0.5609210],
                [-0.4139347, -0.7169559, 0.5609210],
                [-0.5039749, -0.6567929, 0.5609210],
                [-0.5853920, -0.5853920, 0.5609210],
                [-0.6567929, -0.5039749, 0.5609210],
                [-0.7169559, -0.4139347, 0.5609210],
                [-0.7648516, -0.3168119, 0.5609210],
                [-0.7996604, -0.2142684, 0.5609210],
                [-0.8207868, -0.1080586, 0.5609210],
                [-0.8278694, 0.0000000, 0.5609210],
                [-0.8207868, 0.1080586, 0.5609210],
                [-0.7996604, 0.2142684, 0.5609210],
                [-0.7648516, 0.3168119, 0.5609210],
                [-0.7169559, 0.4139347, 0.5609210],
                [-0.6567929, 0.5039749, 0.5609210],
                [-0.5853920, 0.5853920, 0.5609210],
                [-0.5039749, 0.6567929, 0.5609210],
                [-0.4139347, 0.7169559, 0.5609210],
                [-0.3168119, 0.7648516, 0.5609210],
                [-0.2142684, 0.7996604, 0.5609210],
                [-0.1080586, 0.8207868, 0.5609210],
                [0.0000000, 0.7624648, 0.6470297],
                [0.0995216, 0.7559418, 0.6470297],
                [0.1973404, 0.7364845, 0.6470297],
                [0.2917826, 0.7044256, 0.6470297],
                [0.3812324, 0.6603139, 0.6470297],
                [0.4641592, 0.6049040, 0.6470297],
                [0.5391440, 0.5391440, 0.6470297],
                [0.6049040, 0.4641592, 0.6470297],
                [0.6603139, 0.3812324, 0.6470297],
                [0.7044256, 0.2917826, 0.6470297],
                [0.7364845, 0.1973404, 0.6470297],
                [0.7559418, 0.0995216, 0.6470297],
                [0.7624648, 0.0000000, 0.6470297],
                [0.7559418, -0.0995216, 0.6470297],
                [0.7364845, -0.1973404, 0.6470297],
                [0.7044256, -0.2917826, 0.6470297],
                [0.6603139, -0.3812324, 0.6470297],
                [0.6049040, -0.4641592, 0.6470297],
                [0.5391440, -0.5391440, 0.6470297],
                [0.4641592, -0.6049040, 0.6470297],
                [0.3812324, -0.6603139, 0.6470297],
                [0.2917826, -0.7044256, 0.6470297],
                [0.1973404, -0.7364845, 0.6470297],
                [0.0995216, -0.7559418, 0.6470297],
                [0.0000000, -0.7624648, 0.6470297],
                [-0.0995216, -0.7559418, 0.6470297],
                [-0.1973404, -0.7364845, 0.6470297],
                [-0.2917826, -0.7044256, 0.6470297],
                [-0.3812324, -0.6603139, 0.6470297],
                [-0.4641592, -0.6049040, 0.6470297],
                [-0.5391440, -0.5391440, 0.6470297],
                [-0.6049040, -0.4641592, 0.6470297],
                [-0.6603139, -0.3812324, 0.6470297],
                [-0.7044256, -0.2917826, 0.6470297],
                [-0.7364845, -0.1973404, 0.6470297],
                [-0.7559418, -0.0995216, 0.6470297],
                [-0.7624648, 0.0000000, 0.6470297],
                [-0.7559418, 0.0995216, 0.6470297],
                [-0.7364845, 0.1973404, 0.6470297],
                [-0.7044256, 0.2917826, 0.6470297],
                [-0.6603139, 0.3812324, 0.6470297],
                [-0.6049040, 0.4641592, 0.6470297],
                [-0.5391440, 0.5391440, 0.6470297],
                [-0.4641592, 0.6049040, 0.6470297],
                [-0.3812324, 0.6603139, 0.6470297],
                [-0.2917826, 0.7044256, 0.6470297],
                [-0.1973404, 0.7364845, 0.6470297],
                [-0.0995216, 0.7559418, 0.6470297],
                [0.0000000, 0.6881822, 0.7255379],
                [0.0898258, 0.6822947, 0.7255379],
                [0.1781147, 0.6647330, 0.7255379],
                [0.2633559, 0.6357975, 0.7255379],
                [0.3440911, 0.5959833, 0.7255379],
                [0.4189388, 0.5459717, 0.7255379],
                [0.4866183, 0.4866183, 0.7255379],
                [0.5459717, 0.4189388, 0.7255379],
                [0.5959833, 0.3440911, 0.7255379],
                [0.6357975, 0.2633559, 0.7255379],
                [0.6647330, 0.1781147, 0.7255379],
                [0.6822947, 0.0898258, 0.7255379],
                [0.6881822, 0.0000000, 0.7255379],
                [0.6822947, -0.0898258, 0.7255379],
                [0.6647330, -0.1781147, 0.7255379],
                [0.6357975, -0.2633559, 0.7255379],
                [0.5959833, -0.3440911, 0.7255379],
                [0.5459717, -0.4189388, 0.7255379],
                [0.4866183, -0.4866183, 0.7255379],
                [0.4189388, -0.5459717, 0.7255379],
                [0.3440911, -0.5959833, 0.7255379],
                [0.2633559, -0.6357975, 0.7255379],
                [0.1781147, -0.6647330, 0.7255379],
                [0.0898258, -0.6822947, 0.7255379],
                [0.0000000, -0.6881822, 0.7255379],
                [-0.0898258, -0.6822947, 0.7255379],
                [-0.1781147, -0.6647330, 0.7255379],
                [-0.2633559, -0.6357975, 0.7255379],
                [-0.3440911, -0.5959833, 0.7255379],
                [-0.4189388, -0.5459717, 0.7255379],
                [-0.4866183, -0.4866183, 0.7255379],
                [-0.5459717, -0.4189388, 0.7255379],
                [-0.5959833, -0.3440911, 0.7255379],
                [-0.6357975, -0.2633559, 0.7255379],
                [-0.6647330, -0.1781147, 0.7255379],
                [-0.6822947, -0.0898258, 0.7255379],
                [-0.6881822, 0.0000000, 0.7255379],
                [-0.6822947, 0.0898258, 0.7255379],
                [-0.6647330, 0.1781147, 0.7255379],
                [-0.6357975, 0.2633559, 0.7255379],
                [-0.5959833, 0.3440911, 0.7255379],
                [-0.5459717, 0.4189388, 0.7255379],
                [-0.4866183, 0.4866183, 0.7255379],
                [-0.4189388, 0.5459717, 0.7255379],
                [-0.3440911, 0.5959833, 0.7255379],
                [-0.2633559, 0.6357975, 0.7255379],
                [-0.1781147, 0.6647330, 0.7255379],
                [-0.0898258, 0.6822947, 0.7255379],
                [0.0000000, 0.6056963, 0.7956959],
                [0.1051781, 0.5964944, 0.7956959],
                [0.2071603, 0.5691683, 0.7956959],
                [0.3028481, 0.5245484, 0.7956959],
                [0.3893341, 0.4639903, 0.7956959],
                [0.4639903, 0.3893341, 0.7956959],
                [0.5245484, 0.3028481, 0.7956959],
                [0.5691683, 0.2071603, 0.7956959],
                [0.5964944, 0.1051781, 0.7956959],
                [0.6056963, 0.0000000, 0.7956959],
                [0.5964944, -0.1051781, 0.7956959],
                [0.5691683, -0.2071603, 0.7956959],
                [0.5245484, -0.3028481, 0.7956959],
                [0.4639903, -0.3893341, 0.7956959],
                [0.3893341, -0.4639903, 0.7956959],
                [0.3028481, -0.5245484, 0.7956959],
                [0.2071603, -0.5691683, 0.7956959],
                [0.1051781, -0.5964944, 0.7956959],
                [0.0000000, -0.6056963, 0.7956959],
                [-0.1051781, -0.5964944, 0.7956959],
                [-0.2071603, -0.5691683, 0.7956959],
                [-0.3028481, -0.5245484, 0.7956959],
                [-0.3893341, -0.4639903, 0.7956959],
                [-0.4639903, -0.3893341, 0.7956959],
                [-0.5245484, -0.3028481, 0.7956959],
                [-0.5691683, -0.2071603, 0.7956959],
                [-0.5964944, -0.1051781, 0.7956959],
                [-0.6056963, 0.0000000, 0.7956959],
                [-0.5964944, 0.1051781, 0.7956959],
                [-0.5691683, 0.2071603, 0.7956959],
                [-0.5245484, 0.3028481, 0.7956959],
                [-0.4639903, 0.3893341, 0.7956959],
                [-0.3893341, 0.4639903, 0.7956959],
                [-0.3028481, 0.5245484, 0.7956959],
                [-0.2071603, 0.5691683, 0.7956959],
                [-0.1051781, 0.5964944, 0.7956959],
                [0.0000000, 0.5164506, 0.8563170],
                [0.0896807, 0.5086046, 0.8563170],
                [0.1766365, 0.4853048, 0.8563170],
                [0.2582253, 0.4472594, 0.8563170],
                [0.3319681, 0.3956241, 0.8563170],
                [0.3956241, 0.3319681, 0.8563170],
                [0.4472594, 0.2582253, 0.8563170],
                [0.4853048, 0.1766365, 0.8563170],
                [0.5086046, 0.0896807, 0.8563170],
                [0.5164506, 0.0000000, 0.8563170],
                [0.5086046, -0.0896807, 0.8563170],
                [0.4853048, -0.1766365, 0.8563170],
                [0.4472594, -0.2582253, 0.8563170],
                [0.3956241, -0.3319681, 0.8563170],
                [0.3319681, -0.3956241, 0.8563170],
                [0.2582253, -0.4472594, 0.8563170],
                [0.1766365, -0.4853048, 0.8563170],
                [0.0896807, -0.5086046, 0.8563170],
                [0.0000000, -0.5164506, 0.8563170],
                [-0.0896807, -0.5086046, 0.8563170],
                [-0.1766365, -0.4853048, 0.8563170],
                [-0.2582253, -0.4472594, 0.8563170],
                [-0.3319681, -0.3956241, 0.8563170],
                [-0.3956241, -0.3319681, 0.8563170],
                [-0.4472594, -0.2582253, 0.8563170],
                [-0.4853048, -0.1766365, 0.8563170],
                [-0.5086046, -0.0896807, 0.8563170],
                [-0.5164506, 0.0000000, 0.8563170],
                [-0.5086046, 0.0896807, 0.8563170],
                [-0.4853048, 0.1766365, 0.8563170],
                [-0.4472594, 0.2582253, 0.8563170],
                [-0.3956241, 0.3319681, 0.8563170],
                [-0.3319681, 0.3956241, 0.8563170],
                [-0.2582253, 0.4472594, 0.8563170],
                [-0.1766365, 0.4853048, 0.8563170],
                [-0.0896807, 0.5086046, 0.8563170],
                [0.0000000, 0.4208095, 0.9071490],
                [0.1089135, 0.4064708, 0.9071490],
                [0.2104047, 0.3644317, 0.9071490],
                [0.2975572, 0.2975572, 0.9071490],
                [0.3644317, 0.2104047, 0.9071490],
                [0.4064708, 0.1089135, 0.9071490],
                [0.4208095, 0.0000000, 0.9071490],
                [0.4064708, -0.1089135, 0.9071490],
                [0.3644317, -0.2104047, 0.9071490],
                [0.2975572, -0.2975572, 0.9071490],
                [0.2104047, -0.3644317, 0.9071490],
                [0.1089135, -0.4064708, 0.9071490],
                [0.0000000, -0.4208095, 0.9071490],
                [-0.1089135, -0.4064708, 0.9071490],
                [-0.2104047, -0.3644317, 0.9071490],
                [-0.2975572, -0.2975572, 0.9071490],
                [-0.3644317, -0.2104047, 0.9071490],
                [-0.4064708, -0.1089135, 0.9071490],
                [-0.4208095, 0.0000000, 0.9071490],
                [-0.4064708, 0.1089135, 0.9071490],
                [-0.3644317, 0.2104047, 0.9071490],
                [-0.2975572, 0.2975572, 0.9071490],
                [-0.2104047, 0.3644317, 0.9071490],
                [-0.1089135, 0.4064708, 0.9071490],
                [0.0000000, 0.3212214, 0.9470041],
                [0.0831382, 0.3102761, 0.9470041],
                [0.1606107, 0.2781859, 0.9470041],
                [0.2271378, 0.2271378, 0.9470041],
                [0.2781859, 0.1606107, 0.9470041],
                [0.3102761, 0.0831382, 0.9470041],
                [0.3212214, 0.0000000, 0.9470041],
                [0.3102761, -0.0831382, 0.9470041],
                [0.2781859, -0.1606107, 0.9470041],
                [0.2271378, -0.2271378, 0.9470041],
                [0.1606107, -0.2781859, 0.9470041],
                [0.0831382, -0.3102761, 0.9470041],
                [0.0000000, -0.3212214, 0.9470041],
                [-0.0831382, -0.3102761, 0.9470041],
                [-0.1606107, -0.2781859, 0.9470041],
                [-0.2271378, -0.2271378, 0.9470041],
                [-0.2781859, -0.1606107, 0.9470041],
                [-0.3102761, -0.0831382, 0.9470041],
                [-0.3212214, 0.0000000, 0.9470041],
                [-0.3102761, 0.0831382, 0.9470041],
                [-0.2781859, 0.1606107, 0.9470041],
                [-0.2271378, 0.2271378, 0.9470041],
                [-0.1606107, 0.2781859, 0.9470041],
                [-0.0831382, 0.3102761, 0.9470041],
                [0.0000000, 0.2169121, 0.9761911],
                [0.1084561, 0.1878514, 0.9761911],
                [0.1878514, 0.1084561, 0.9761911],
                [0.2169121, 0.0000000, 0.9761911],
                [0.1878514, -0.1084561, 0.9761911],
                [0.1084561, -0.1878514, 0.9761911],
                [0.0000000, -0.2169121, 0.9761911],
                [-0.1084561, -0.1878514, 0.9761911],
                [-0.1878514, -0.1084561, 0.9761911],
                [-0.2169121, 0.0000000, 0.9761911],
                [-0.1878514, 0.1084561, 0.9761911],
                [-0.1084561, 0.1878514, 0.9761911],
                [0.0000000, 0.1157325, 0.9932804],
                [0.0578662, 0.1002273, 0.9932804],
                [0.1002273, 0.0578662, 0.9932804],
                [0.1157325, 0.0000000, 0.9932804],
                [0.1002273, -0.0578662, 0.9932804],
                [0.0578662, -0.1002273, 0.9932804],
                [0.0000000, -0.1157325, 0.9932804],
                [-0.0578662, -0.1002273, 0.9932804],
                [-0.1002273, -0.0578662, 0.9932804],
                [-0.1157325, 0.0000000, 0.9932804],
                [-0.1002273, 0.0578662, 0.9932804],
                [-0.0578662, 0.1002273, 0.9932804],
                [0.0000000, -0.0000000, 1.0000000]
            ])
            self.patch_count = 577
            self.row_conversion_factor = np.array([
                0.0113221971,
                0.0111894547,
                0.0109255262,
                0.0105335058,
                0.0125224872,
                0.0117312774,
                0.0108025291,
                0.00974713106,
                0.011436609,
                0.00974295956,
                0.0119026242,
                0.00905126163,
                0.0121875626,
                0.00612971396,
                0.00921483254,
            ])
            self.row_patches = np.array([
                60,
                60,
                60,
                60,
                48,
                48,
                48,
                48,
                36,
                36,
                24,
                24,
                12,
                12,
                1,
            ])
            self.patch_conversion_factor = np.array([
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.0113222, 0.0113222, 0.0113222, 0.0113222, 0.0113222,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01118945, 0.01118945, 0.01118945, 0.01118945, 0.01118945,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01092553, 0.01092553, 0.01092553, 0.01092553, 0.01092553,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01053351, 0.01053351, 0.01053351, 0.01053351, 0.01053351,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01252249, 0.01252249,
                0.01252249, 0.01252249, 0.01252249, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01173128, 0.01173128, 0.01173128, 0.01173128,
                0.01173128, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.01080253,
                0.01080253, 0.01080253, 0.01080253, 0.01080253, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.00974713, 0.00974713, 0.00974713,
                0.00974713, 0.00974713, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.01143661, 0.01143661,
                0.01143661, 0.01143661, 0.01143661, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.00974296,
                0.00974296, 0.00974296, 0.00974296, 0.00974296, 0.01190262,
                0.01190262, 0.01190262, 0.01190262, 0.01190262, 0.01190262,
                0.01190262, 0.01190262, 0.01190262, 0.01190262, 0.01190262,
                0.01190262, 0.01190262, 0.01190262, 0.01190262, 0.01190262,
                0.01190262, 0.01190262, 0.01190262, 0.01190262, 0.01190262,
                0.01190262, 0.01190262, 0.01190262, 0.00905126, 0.00905126,
                0.00905126, 0.00905126, 0.00905126, 0.00905126, 0.00905126,
                0.00905126, 0.00905126, 0.00905126, 0.00905126, 0.00905126,
                0.00905126, 0.00905126, 0.00905126, 0.00905126, 0.00905126,
                0.00905126, 0.00905126, 0.00905126, 0.00905126, 0.00905126,
                0.00905126, 0.00905126, 0.01218756, 0.01218756, 0.01218756,
                0.01218756, 0.01218756, 0.01218756, 0.01218756, 0.01218756,
                0.01218756, 0.01218756, 0.01218756, 0.01218756, 0.00612971,
                0.00612971, 0.00612971, 0.00612971, 0.00612971, 0.00612971,
                0.00612971, 0.00612971, 0.00612971, 0.00612971, 0.00612971,
                0.00612971, 0.00921483
            ])
        else:
            self.patch_centroids = np.array([
                [0, 99.756283, 5.405868],
                [10.427371, 99.209807, 5.405868],
                [20.740497, 97.576367, 5.405868],
                [30.826387, 94.873862, 5.405868],
                [40.574535, 91.131899, 5.405868],
                [49.878141, 86.391475, 5.405868],
                [58.635271, 80.704528, 5.405868],
                [66.749981, 74.133365, 5.405868],
                [74.133365, 66.749981, 5.405868],
                [80.704528, 58.635271, 5.405868],
                [86.391475, 49.878141, 5.405868],
                [91.131899, 40.574535, 5.405868],
                [94.873862, 30.826387, 5.405868],
                [97.576367, 20.740497, 5.405868],
                [99.209807, 10.427371, 5.405868],
                [99.756283, 0, 5.405868],
                [99.209807, -10.427371, 5.405868],
                [97.576367, -20.740497, 5.405868],
                [94.873862, -30.826387, 5.405868],
                [91.131899, -40.574535, 5.405868],
                [86.391475, -49.878141, 5.405868],
                [80.704528, -58.635271, 5.405868],
                [74.133365, -66.749981, 5.405868],
                [66.749981, -74.133365, 5.405868],
                [58.635271, -80.704528, 5.405868],
                [49.878141, -86.391475, 5.405868],
                [40.574535, -91.131899, 5.405868],
                [30.826387, -94.873862, 5.405868],
                [20.740497, -97.576367, 5.405868],
                [10.427371, -99.209807, 5.405868],
                [0, -99.756283, 5.405868],
                [-10.427371, -99.209807, 5.405868],
                [-20.740497, -97.576367, 5.405868],
                [-30.826387, -94.873862, 5.405868],
                [-40.574535, -91.131899, 5.405868],
                [-49.878141, -86.391475, 5.405868],
                [-58.635271, -80.704528, 5.405868],
                [-66.749981, -74.133365, 5.405868],
                [-74.133365, -66.749981, 5.405868],
                [-80.704528, -58.635271, 5.405868],
                [-86.391475, -49.878141, 5.405868],
                [-91.131899, -40.574535, 5.405868],
                [-94.873862, -30.826387, 5.405868],
                [-97.576367, -20.740497, 5.405868],
                [-99.209807, -10.427371, 5.405868],
                [-99.756283, 0, 5.405868],
                [-99.209807, 10.427371, 5.405868],
                [-97.576367, 20.740497, 5.405868],
                [-94.873862, 30.826387, 5.405868],
                [-91.131899, 40.574535, 5.405868],
                [-86.391475, 49.878141, 5.405868],
                [-80.704528, 58.635271, 5.405868],
                [-74.133365, 66.749981, 5.405868],
                [-66.749981, 74.133365, 5.405868],
                [-58.635271, 80.704528, 5.405868],
                [-49.878141, 86.391475, 5.405868],
                [-40.574535, 91.131899, 5.405868],
                [-30.826387, 94.873862, 5.405868],
                [-20.740497, 97.576367, 5.405868],
                [-10.427371, 99.209807, 5.405868],
                [0, 98.58904, 16.154225],
                [10.305361, 98.048959, 16.154225],
                [20.497814, 96.434632, 16.154225],
                [30.465689, 93.763749, 16.154226],
                [40.099775, 90.06557, 16.154225],
                [49.29452, 85.380613, 16.154225],
                [57.949183, 79.760209, 16.154226],
                [65.968943, 73.265935, 16.154225],
                [73.265935, 65.968943, 16.154225],
                [79.760209, 57.949183, 16.154226],
                [85.380613, 49.29452, 16.154225],
                [90.06557, 40.099775, 16.154225],
                [93.763749, 30.465689, 16.154226],
                [96.434632, 20.497814, 16.154225],
                [98.048959, 10.305361, 16.154225],
                [98.58904, 0, 16.154225],
                [98.048959, -10.305361, 16.154225],
                [96.434632, -20.497814, 16.154225],
                [93.763749, -30.465689, 16.154226],
                [90.06557, -40.099775, 16.154225],
                [85.380613, -49.29452, 16.154225],
                [79.760209, -57.949183, 16.154226],
                [73.265935, -65.968943, 16.154225],
                [65.968943, -73.265935, 16.154225],
                [57.949183, -79.760209, 16.154226],
                [49.29452, -85.380613, 16.154225],
                [40.099775, -90.06557, 16.154225],
                [30.465689, -93.763749, 16.154226],
                [20.497814, -96.434632, 16.154225],
                [10.305361, -98.048959, 16.154225],
                [0, -98.58904, 16.154225],
                [-10.305361, -98.048959, 16.154225],
                [-20.497814, -96.434632, 16.154225],
                [-30.465689, -93.763749, 16.154226],
                [-40.099775, -90.06557, 16.154225],
                [-49.29452, -85.380613, 16.154225],
                [-57.949183, -79.760209, 16.154226],
                [-65.968943, -73.265935, 16.154225],
                [-73.265935, -65.968943, 16.154225],
                [-79.760209, -57.949183, 16.154226],
                [-85.380613, -49.29452, 16.154225],
                [-90.06557, -40.099775, 16.154225],
                [-93.763749, -30.465689, 16.154226],
                [-96.434632, -20.497814, 16.154225],
                [-98.048959, -10.305361, 16.154225],
                [-98.58904, 0, 16.154225],
                [-98.048959, 10.305361, 16.154225],
                [-96.434632, 20.497814, 16.154225],
                [-93.763749, 30.465689, 16.154226],
                [-90.06557, 40.099775, 16.154225],
                [-85.380613, 49.29452, 16.154225],
                [-79.760209, 57.949183, 16.154226],
                [-73.265935, 65.968943, 16.154225],
                [-65.968943, 73.265935, 16.154225],
                [-57.949183, 79.760209, 16.154226],
                [-49.29452, 85.380613, 16.154225],
                [-40.099775, 90.06557, 16.154225],
                [-30.465689, 93.763749, 16.154226],
                [-20.497814, 96.434632, 16.154225],
                [-10.305361, 98.048959, 16.154225],
                [0, 96.268323, 26.713189],
                [10.06278, 95.740956, 26.713189],
                [20.01531, 94.164631, 26.713189],
                [29.748548, 91.556616, 26.713189],
                [39.155855, 87.94549, 26.713189],
                [48.134162, 83.370815, 26.713189],
                [56.585101, 77.88271, 26.713189],
                [64.416082, 71.541307, 26.713189],
                [71.541307, 64.416082, 26.713189],
                [77.88271, 56.585101, 26.713189],
                [83.370815, 48.134162, 26.713189],
                [87.94549, 39.155855, 26.713189],
                [91.556616, 29.748548, 26.713189],
                [94.164631, 20.01531, 26.713189],
                [95.740956, 10.06278, 26.713189],
                [96.268323, 0, 26.713189],
                [95.740956, -10.06278, 26.713189],
                [94.164631, -20.01531, 26.713189],
                [91.556616, -29.748548, 26.713189],
                [87.94549, -39.155855, 26.713189],
                [83.370815, -48.134162, 26.713189],
                [77.88271, -56.585101, 26.713189],
                [71.541307, -64.416082, 26.713189],
                [64.416082, -71.541307, 26.713189],
                [56.585101, -77.88271, 26.713189],
                [48.134162, -83.370815, 26.713189],
                [39.155855, -87.94549, 26.713189],
                [29.748548, -91.556616, 26.713189],
                [20.01531, -94.164631, 26.713189],
                [10.06278, -95.740956, 26.713189],
                [0, -96.268323, 26.713189],
                [-10.06278, -95.740956, 26.713189],
                [-20.01531, -94.164631, 26.713189],
                [-29.748548, -91.556616, 26.713189],
                [-39.155855, -87.94549, 26.713189],
                [-48.134162, -83.370815, 26.713189],
                [-56.585101, -77.88271, 26.713189],
                [-64.416082, -71.541307, 26.713189],
                [-71.541307, -64.416082, 26.713189],
                [-77.88271, -56.585101, 26.713189],
                [-83.370815, -48.134162, 26.713189],
                [-87.94549, -39.155855, 26.713189],
                [-91.556616, -29.748548, 26.713189],
                [-94.164631, -20.01531, 26.713189],
                [-95.740956, -10.06278, 26.713189],
                [-96.268323, 0, 26.713189],
                [-95.740956, 10.06278, 26.713189],
                [-94.164631, 20.01531, 26.713189],
                [-91.556616, 29.748548, 26.713189],
                [-87.94549, 39.155855, 26.713189],
                [-83.370815, 48.134162, 26.713189],
                [-77.88271, 56.585101, 26.713189],
                [-71.541307, 64.416082, 26.713189],
                [-64.416082, 71.541307, 26.713189],
                [-56.585101, 77.88271, 26.713189],
                [-48.134162, 83.370815, 26.713189],
                [-39.155855, 87.94549, 26.713189],
                [-29.748548, 91.556616, 26.713189],
                [-20.01531, 94.164631, 26.713189],
                [-10.06278, 95.740956, 26.713189],
                [0, 92.821522, 36.958966],
                [9.702491, 92.313036, 36.958966],
                [19.29868, 90.793149, 36.958966],
                [28.683428, 88.278514, 36.958966],
                [37.753914, 84.79668, 36.958966],
                [46.410761, 80.385796, 36.958966],
                [54.559121, 75.094189, 36.958966],
                [62.109721, 68.979834, 36.958966],
                [68.979834, 62.109721, 36.958966],
                [75.094189, 54.559121, 36.958966],
                [80.385796, 46.410761, 36.958966],
                [84.79668, 37.753914, 36.958966],
                [88.278514, 28.683428, 36.958966],
                [90.793149, 19.29868, 36.958966],
                [92.313036, 9.702491, 36.958966],
                [92.821522, 0, 36.958966],
                [92.313036, -9.702491, 36.958966],
                [90.793149, -19.29868, 36.958966],
                [88.278514, -28.683428, 36.958966],
                [84.79668, -37.753914, 36.958966],
                [80.385796, -46.410761, 36.958966],
                [75.094189, -54.559121, 36.958966],
                [68.979834, -62.109721, 36.958966],
                [62.109721, -68.979834, 36.958966],
                [54.559121, -75.094189, 36.958966],
                [46.410761, -80.385796, 36.958966],
                [37.753914, -84.79668, 36.958966],
                [28.683428, -88.278514, 36.958966],
                [19.29868, -90.793149, 36.958966],
                [9.702491, -92.313036, 36.958966],
                [0, -92.821522, 36.958966],
                [-9.702491, -92.313036, 36.958966],
                [-19.29868, -90.793149, 36.958966],
                [-28.683428, -88.278514, 36.958966],
                [-37.753914, -84.79668, 36.958966],
                [-46.410761, -80.385796, 36.958966],
                [-54.559121, -75.094189, 36.958966],
                [-62.109721, -68.979834, 36.958966],
                [-68.979834, -62.109721, 36.958966],
                [-75.094189, -54.559121, 36.958966],
                [-80.385796, -46.410761, 36.958966],
                [-84.79668, -37.753914, 36.958966],
                [-88.278514, -28.683428, 36.958966],
                [-90.793149, -19.29868, 36.958966],
                [-92.313036, -9.702491, 36.958966],
                [-92.821522, 0, 36.958966],
                [-92.313036, 9.702491, 36.958966],
                [-90.793149, 19.29868, 36.958966],
                [-88.278514, 28.683428, 36.958966],
                [-84.79668, 37.753914, 36.958966],
                [-80.385796, 46.410761, 36.958966],
                [-75.094189, 54.559121, 36.958966],
                [-68.979834, 62.109721, 36.958966],
                [-62.109721, 68.979834, 36.958966],
                [-54.559121, 75.094189, 36.958966],
                [-46.410761, 80.385796, 36.958966],
                [-37.753914, 84.79668, 36.958966],
                [-28.683428, 88.278514, 36.958966],
                [-19.29868, 90.793149, 36.958966],
                [-9.702491, 92.313036, 36.958966],
                [0, 88.265941, 46.771431],
                [11.521017, 87.510814, 46.771431],
                [22.844907, 85.258352, 46.771431],
                [33.777914, 81.547096, 46.771431],
                [44.132971, 76.440547, 46.771431],
                [53.732901, 70.026079, 46.771431],
                [62.413445, 62.413445, 46.771431],
                [70.026079, 53.732901, 46.771431],
                [76.440547, 44.132971, 46.771431],
                [81.547096, 33.777914, 46.771431],
                [85.258352, 22.844907, 46.771431],
                [87.510814, 11.521017, 46.771431],
                [88.265941, 0, 46.771431],
                [87.510814, -11.521017, 46.771431],
                [85.258352, -22.844907, 46.771431],
                [81.547096, -33.777914, 46.771431],
                [76.440547, -44.132971, 46.771431],
                [70.026079, -53.732901, 46.771431],
                [62.413445, -62.413445, 46.771431],
                [53.732901, -70.026079, 46.771431],
                [44.132971, -76.440547, 46.771431],
                [33.777914, -81.547096, 46.771431],
                [22.844907, -85.258352, 46.771431],
                [11.521017, -87.510814, 46.771431],
                [0, -88.265941, 46.771431],
                [-11.521017, -87.510814, 46.771431],
                [-22.844907, -85.258352, 46.771431],
                [-33.777914, -81.547096, 46.771431],
                [-44.132971, -76.440547, 46.771431],
                [-53.732901, -70.026079, 46.771431],
                [-62.413445, -62.413445, 46.771431],
                [-70.026079, -53.732901, 46.771431],
                [-76.440547, -44.132971, 46.771431],
                [-81.547096, -33.777914, 46.771431],
                [-85.258352, -22.844907, 46.771431],
                [-87.510814, -11.521017, 46.771431],
                [-88.265941, 0, 46.771431],
                [-87.510814, 11.521017, 46.771431],
                [-85.258352, 22.844907, 46.771431],
                [-81.547096, 33.777914, 46.771431],
                [-76.440547, 44.132971, 46.771431],
                [-70.026079, 53.732901, 46.771431],
                [-62.413445, 62.413445, 46.771431],
                [-53.732901, 70.026079, 46.771431],
                [-44.132971, 76.440547, 46.771431],
                [-33.777914, 81.547096, 46.771431],
                [-22.844907, 85.258352, 46.771431],
                [-11.521017, 87.510814, 46.771431],
                [0, 82.703471, 56.035545],
                [10.794969, 81.995932, 56.035545],
                [21.405233, 79.885419, 56.035545],
                [31.649248, 76.408044, 56.035545],
                [41.351736, 71.623307, 56.035545],
                [50.346684, 65.613075, 56.035545],
                [58.480186, 58.480186, 56.035545],
                [65.613075, 50.346684, 56.035545],
                [71.623307, 41.351736, 56.035545],
                [76.408044, 31.649248, 56.035545],
                [79.885419, 21.405233, 56.035545],
                [81.995932, 10.794969, 56.035545],
                [82.703471, 0, 56.035545],
                [81.995932, -10.794969, 56.035545],
                [79.885419, -21.405233, 56.035545],
                [76.408044, -31.649248, 56.035545],
                [71.623307, -41.351736, 56.035545],
                [65.613075, -50.346684, 56.035545],
                [58.480186, -58.480186, 56.035545],
                [50.346684, -65.613075, 56.035545],
                [41.351736, -71.623307, 56.035545],
                [31.649248, -76.408044, 56.035545],
                [21.405233, -79.885419, 56.035545],
                [10.794969, -81.995932, 56.035545],
                [0, -82.703471, 56.035545],
                [-10.794969, -81.995932, 56.035545],
                [-21.405233, -79.885419, 56.035545],
                [-31.649248, -76.408044, 56.035545],
                [-41.351736, -71.623307, 56.035545],
                [-50.346684, -65.613075, 56.035545],
                [-58.480186, -58.480186, 56.035545],
                [-65.613075, -50.346684, 56.035545],
                [-71.623307, -41.351736, 56.035545],
                [-76.408044, -31.649248, 56.035545],
                [-79.885419, -21.405233, 56.035545],
                [-81.995932, -10.794969, 56.035545],
                [-82.703471, 0, 56.035545],
                [-81.995932, 10.794969, 56.035545],
                [-79.885419, 21.405233, 56.035545],
                [-76.408044, 31.649248, 56.035545],
                [-71.623307, 41.351736, 56.035545],
                [-65.613075, 50.346684, 56.035545],
                [-58.480186, 58.480186, 56.035545],
                [-50.346684, 65.613075, 56.035545],
                [-41.351736, 71.623307, 56.035545],
                [-31.649248, 76.408044, 56.035545],
                [-21.405233, 79.885419, 56.035545],
                [-10.794969, 81.995932, 56.035545],
                [0, 76.175454, 64.642694],
                [9.942892, 75.523762, 64.642693],
                [19.715658, 73.579839, 64.642693],
                [29.151084, 70.376943, 64.642693],
                [38.087727, 65.969879, 64.642693],
                [46.372678, 60.434051, 64.642693],
                [53.864181, 53.864181, 64.642692],
                [60.434051, 46.372678, 64.642693],
                [65.969879, 38.087727, 64.642693],
                [70.376943, 29.151084, 64.642693],
                [73.579839, 19.715658, 64.642693],
                [75.523762, 9.942892, 64.642693],
                [76.175454, 0, 64.642694],
                [75.523762, -9.942892, 64.642693],
                [73.579839, -19.715658, 64.642693],
                [70.376943, -29.151084, 64.642693],
                [65.969879, -38.087727, 64.642693],
                [60.434051, -46.372678, 64.642693],
                [53.864181, -53.864181, 64.642692],
                [46.372678, -60.434051, 64.642693],
                [38.087727, -65.969879, 64.642693],
                [29.151084, -70.376943, 64.642693],
                [19.715658, -73.579839, 64.642693],
                [9.942892, -75.523762, 64.642693],
                [0, -76.175454, 64.642694],
                [-9.942892, -75.523762, 64.642693],
                [-19.715658, -73.579839, 64.642693],
                [-29.151084, -70.376943, 64.642693],
                [-38.087727, -65.969879, 64.642693],
                [-46.372678, -60.434051, 64.642693],
                [-53.864181, -53.864181, 64.642692],
                [-60.434051, -46.372678, 64.642693],
                [-65.969879, -38.087727, 64.642693],
                [-70.376943, -29.151084, 64.642693],
                [-73.579839, -19.715658, 64.642693],
                [-75.523762, -9.942892, 64.642693],
                [-76.175454, 0, 64.642694],
                [-75.523762, 9.942892, 64.642693],
                [-73.579839, 19.715658, 64.642693],
                [-70.376943, 29.151084, 64.642693],
                [-65.969879, 38.087727, 64.642693],
                [-60.434051, 46.372678, 64.642693],
                [-53.864181, 53.864181, 64.642692],
                [-46.372678, 60.434051, 64.642693],
                [-38.087727, 65.969879, 64.642693],
                [-29.151084, 70.376943, 64.642693],
                [-19.715658, 73.579839, 64.642693],
                [-9.942892, 75.523762, 64.642693],
                [0, 68.759585, 72.491965],
                [8.974927, 68.171338, 72.491965],
                [17.796291, 66.416659, 72.491965],
                [26.313154, 63.525574, 72.491965],
                [34.379793, 59.547548, 72.491965],
                [41.858183, 54.550646, 72.491965],
                [48.62037, 48.62037, 72.491965],
                [54.550646, 41.858183, 72.491965],
                [59.547548, 34.379793, 72.491965],
                [63.525574, 26.313154, 72.491965],
                [66.416659, 17.796291, 72.491965],
                [68.171338, 8.974927, 72.491965],
                [68.759585, 0, 72.491965],
                [68.171338, -8.974927, 72.491965],
                [66.416659, -17.796291, 72.491965],
                [63.525574, -26.313154, 72.491965],
                [59.547548, -34.379793, 72.491965],
                [54.550646, -41.858183, 72.491965],
                [48.62037, -48.62037, 72.491965],
                [41.858183, -54.550646, 72.491965],
                [34.379793, -59.547548, 72.491965],
                [26.313154, -63.525574, 72.491965],
                [17.796291, -66.416659, 72.491965],
                [8.974927, -68.171338, 72.491965],
                [0, -68.759585, 72.491965],
                [-8.974927, -68.171338, 72.491965],
                [-17.796291, -66.416659, 72.491965],
                [-26.313154, -63.525574, 72.491965],
                [-34.379793, -59.547548, 72.491965],
                [-41.858183, -54.550646, 72.491965],
                [-48.62037, -48.62037, 72.491965],
                [-54.550646, -41.858183, 72.491965],
                [-59.547548, -34.379793, 72.491965],
                [-63.525574, -26.313154, 72.491965],
                [-66.416659, -17.796291, 72.491965],
                [-68.171338, -8.974927, 72.491965],
                [-68.759585, 0, 72.491965],
                [-68.171338, 8.974927, 72.491965],
                [-66.416659, 17.796291, 72.491965],
                [-63.525574, 26.313154, 72.491965],
                [-59.547548, 34.379793, 72.491965],
                [-54.550646, 41.858183, 72.491965],
                [-48.62037, 48.62037, 72.491965],
                [-41.858183, 54.550646, 72.491965],
                [-34.379793, 59.547548, 72.491965],
                [-26.313154, 63.525574, 72.491965],
                [-17.796291, 66.416659, 72.491965],
                [-8.974927, 68.171338, 72.491965],
                [0, 60.510059, 79.491334],
                [10.507461, 59.590775, 79.491333],
                [20.695659, 56.860856, 79.491334],
                [30.255029, 52.403248, 79.491334],
                [38.895116, 46.353394, 79.491334],
                [46.353394, 38.895116, 79.491334],
                [52.403248, 30.255029, 79.491334],
                [56.860856, 20.695659, 79.491334],
                [59.590775, 10.507461, 79.491333],
                [60.510059, 0, 79.491334],
                [59.590775, -10.507461, 79.491333],
                [56.860856, -20.695659, 79.491334],
                [52.403248, -30.255029, 79.491334],
                [46.353394, -38.895116, 79.491334],
                [38.895116, -46.353394, 79.491334],
                [30.255029, -52.403248, 79.491334],
                [20.695659, -56.860856, 79.491334],
                [10.507461, -59.590775, 79.491333],
                [0, -60.510059, 79.491334],
                [-10.507461, -59.590775, 79.491333],
                [-20.695659, -56.860856, 79.491334],
                [-30.255029, -52.403248, 79.491334],
                [-38.895116, -46.353394, 79.491334],
                [-46.353394, -38.895116, 79.491334],
                [-52.403248, -30.255029, 79.491334],
                [-56.860856, -20.695659, 79.491334],
                [-59.590775, -10.507461, 79.491333],
                [-60.510059, 0, 79.491334],
                [-59.590775, 10.507461, 79.491333],
                [-56.860856, 20.695659, 79.491334],
                [-52.403248, 30.255029, 79.491334],
                [-46.353394, 38.895116, 79.491334],
                [-38.895116, 46.353394, 79.491334],
                [-30.255029, 52.403248, 79.491334],
                [-20.695659, 56.860856, 79.491334],
                [-10.507461, 59.590775, 79.491333],
                [0, 51.601062, 85.558742],
                [8.96043, 50.817126, 85.558742],
                [17.648603, 48.489138, 85.558742],
                [25.800531, 44.687831, 85.558742],
                [33.168524, 39.528708, 85.558742],
                [39.528708, 33.168524, 85.558742],
                [44.687831, 25.800531, 85.558742],
                [48.489138, 17.648603, 85.558742],
                [50.817126, 8.96043, 85.558742],
                [51.601062, 0, 85.558742],
                [50.817126, -8.96043, 85.558742],
                [48.489138, -17.648603, 85.558742],
                [44.687831, -25.800531, 85.558742],
                [39.528708, -33.168524, 85.558742],
                [33.168524, -39.528708, 85.558742],
                [25.800531, -44.687831, 85.558742],
                [17.648603, -48.489138, 85.558742],
                [8.96043, -50.817126, 85.558742],
                [0, -51.601062, 85.558742],
                [-8.96043, -50.817126, 85.558742],
                [-17.648603, -48.489138, 85.558742],
                [-25.800531, -44.687831, 85.558742],
                [-33.168524, -39.528708, 85.558742],
                [-39.528708, -33.168524, 85.558742],
                [-44.687831, -25.800531, 85.558742],
                [-48.489138, -17.648603, 85.558742],
                [-50.817126, -8.96043, 85.558742],
                [-51.601062, 0, 85.558742],
                [-50.817126, 8.96043, 85.558742],
                [-48.489138, 17.648603, 85.558742],
                [-44.687831, 25.800531, 85.558742],
                [-39.528708, 33.168524, 85.558742],
                [-33.168524, 39.528708, 85.558742],
                [-25.800531, 44.687831, 85.558742],
                [-17.648603, 48.489138, 85.558742],
                [-8.96043, 50.817126, 85.558742],
                [0, 42.038342, 90.623053],
                [10.880323, 40.605919, 90.623053],
                [21.019171, 36.406271, 90.623053],
                [29.725596, 29.725596, 90.623053],
                [36.406271, 21.019171, 90.623053],
                [40.605919, 10.880323, 90.623053],
                [42.038342, 0, 90.623053],
                [40.605919, -10.880323, 90.623053],
                [36.406271, -21.019171, 90.623053],
                [29.725596, -29.725596, 90.623053],
                [21.019171, -36.406271, 90.623053],
                [10.880323, -40.605919, 90.623053],
                [0, -42.038342, 90.623053],
                [-10.880323, -40.605919, 90.623053],
                [-21.019171, -36.406271, 90.623053],
                [-29.725596, -29.725596, 90.623053],
                [-36.406271, -21.019171, 90.623053],
                [-40.605919, -10.880323, 90.623053],
                [-42.038342, 0, 90.623053],
                [-40.605919, 10.880323, 90.623053],
                [-36.406271, 21.019171, 90.623053],
                [-29.725596, 29.725596, 90.623053],
                [-21.019171, 36.406271, 90.623053],
                [-10.880323, 40.605919, 90.623053],
                [0, 32.096525, 94.624889],
                [8.307192, 31.002862, 94.62489],
                [16.048262, 27.796406, 94.624889],
                [22.69567, 22.69567, 94.62489],
                [27.796406, 16.048262, 94.624889],
                [31.002862, 8.307192, 94.62489],
                [32.096525, 0, 94.624889],
                [31.002862, -8.307192, 94.62489],
                [27.796406, -16.048262, 94.624889],
                [22.69567, -22.69567, 94.62489],
                [16.048262, -27.796406, 94.624889],
                [8.307192, -31.002862, 94.62489],
                [0, -32.096525, 94.624889],
                [-8.307192, -31.002862, 94.62489],
                [-16.048262, -27.796406, 94.624889],
                [-22.69567, -22.69567, 94.62489],
                [-27.796406, -16.048262, 94.624889],
                [-31.002862, -8.307192, 94.62489],
                [-32.096525, 0, 94.624889],
                [-31.002862, 8.307192, 94.62489],
                [-27.796406, 16.048262, 94.624889],
                [-22.69567, 22.69567, 94.62489],
                [-16.048262, 27.796406, 94.624889],
                [-8.307192, 31.002862, 94.62489],
                [0, 21.668599, 97.517335],
                [10.834299, 18.765557, 97.517335],
                [18.765557, 10.834299, 97.517335],
                [21.668599, 0, 97.517335],
                [18.765557, -10.834299, 97.517335],
                [10.834299, -18.765557, 97.517335],
                [0, -21.668599, 97.517335],
                [-10.834299, -18.765557, 97.517335],
                [-18.765557, -10.834299, 97.517335],
                [-21.668599, 0, 97.517335],
                [-18.765557, 10.834299, 97.517335],
                [-10.834299, 18.765557, 97.517335],
                [0, 11.566076, 99.26648],
                [5.783038, 10.016516, 99.26648],
                [10.016516, 5.783038, 99.26648],
                [11.566076, 0, 99.26648],
                [10.016516, -5.783038, 99.26648],
                [5.783038, -10.016516, 99.26648],
                [0, -11.566076, 99.26648],
                [-5.783038, -10.016516, 99.26648],
                [-10.016516, -5.783038, 99.26648],
                [-11.566076, 0, 99.26648],
                [-10.016516, 5.783038, 99.26648],
                [-5.783038, 10.016516, 99.26648],
                [0, 0, 99.925143],
            ])
            self.patch_vectors = np.array([
                [0.0, 0.99452200000000002, 0.104528],
                [0.20677300000000001, 0.97278900000000001, 0.104528],
                [0.40450799999999998, 0.90854100000000004, 0.104528],
                [0.584565, 0.80458499999999999, 0.104528],
                [0.73907400000000001, 0.66546499999999997, 0.104528],
                [0.86128099999999996, 0.49726100000000001, 0.104528],
                [0.94584699999999999, 0.30732399999999999, 0.104528],
                [0.98907400000000001, 0.10395600000000001, 0.104528],
                [0.98907400000000001, -0.10395600000000001, 0.104528],
                [0.94584699999999999, -0.30732399999999999, 0.104528],
                [0.86128099999999996, -0.49726100000000001, 0.104528],
                [0.73907400000000001, -0.66546499999999997, 0.104528],
                [0.584565, -0.80458499999999999, 0.104528],
                [0.40450799999999998, -0.90854100000000004, 0.104528],
                [0.20677300000000001, -0.97278900000000001, 0.104528],
                [0.0, -0.99452200000000002, 0.104528],
                [-0.20677300000000001, -0.97278900000000001, 0.104528],
                [-0.40450799999999998, -0.90854100000000004, 0.104528],
                [-0.584565, -0.80458499999999999, 0.104528],
                [-0.73907400000000001, -0.66546499999999997, 0.104528],
                [-0.86128099999999996, -0.49726100000000001, 0.104528],
                [-0.94584699999999999, -0.30732399999999999, 0.104528],
                [-0.98907400000000001, -0.10395600000000001, 0.104528],
                [-0.98907400000000001, 0.10395600000000001, 0.104528],
                [-0.94584699999999999, 0.30732399999999999, 0.104528],
                [-0.86128099999999996, 0.49726100000000001, 0.104528],
                [-0.73907400000000001, 0.66546499999999997, 0.104528],
                [-0.584565, 0.80458499999999999, 0.104528],
                [-0.40450799999999998, 0.90854100000000004, 0.104528],
                [-0.20677300000000001, 0.97278900000000001, 0.104528],
                [0.0, 0.95105700000000004, 0.30901699999999999],
                [0.197736, 0.93027400000000005, 0.30901699999999999],
                [0.38683000000000001, 0.86883299999999997, 0.30901699999999999],
                [0.55901699999999999, 0.76942100000000002, 0.30901699999999999],
                [0.70677299999999998, 0.63638099999999997, 0.30901699999999999],
                [0.82363900000000001, 0.47552800000000001, 0.30901699999999999],
                [0.90450799999999998, 0.29389300000000002, 0.30901699999999999],
                [0.94584699999999999, 0.099412, 0.30901699999999999],
                [0.94584699999999999, -0.099412, 0.30901699999999999],
                [0.90450799999999998, -0.29389300000000002, 0.30901699999999999],
                [0.82363900000000001, -0.47552800000000001, 0.30901699999999999],
                [0.70677299999999998, -0.63638099999999997, 0.30901699999999999],
                [0.55901699999999999, -0.76942100000000002, 0.30901699999999999],
                [0.38683000000000001, -0.86883299999999997, 0.30901699999999999],
                [0.197736, -0.93027400000000005, 0.30901699999999999],
                [0.0, -0.95105700000000004, 0.30901699999999999],
                [-0.197736, -0.93027400000000005, 0.30901699999999999],
                [-0.38683000000000001, -0.86883299999999997, 0.30901699999999999],
                [-0.55901699999999999, -0.76942100000000002, 0.30901699999999999],
                [-0.70677299999999998, -0.63638099999999997, 0.30901699999999999],
                [-0.82363900000000001, -0.47552800000000001, 0.30901699999999999],
                [-0.90450799999999998, -0.29389300000000002, 0.30901699999999999],
                [-0.94584699999999999, -0.099412, 0.30901699999999999],
                [-0.94584699999999999, 0.099412, 0.30901699999999999],
                [-0.90450799999999998, 0.29389300000000002, 0.30901699999999999],
                [-0.82363900000000001, 0.47552800000000001, 0.30901699999999999],
                [-0.70677299999999998, 0.63638099999999997, 0.30901699999999999],
                [-0.55901699999999999, 0.76942100000000002, 0.30901699999999999],
                [-0.38683000000000001, 0.86883299999999997, 0.30901699999999999],
                [-0.197736, 0.93027400000000005, 0.30901699999999999],
                [0.0, 0.86602500000000004, 0.5],
                [0.22414400000000001, 0.83651600000000004, 0.5],
                [0.43301299999999998, 0.75, 0.5],
                [0.61237200000000003, 0.61237200000000003, 0.5],
                [0.75, 0.43301299999999998, 0.5],
                [0.83651600000000004, 0.22414400000000001, 0.5],
                [0.86602500000000004, 0.0, 0.5],
                [0.83651600000000004, -0.22414400000000001, 0.5],
                [0.75, -0.43301299999999998, 0.5],
                [0.61237200000000003, -0.61237200000000003, 0.5],
                [0.43301299999999998, -0.75, 0.5],
                [0.22414400000000001, -0.83651600000000004, 0.5],
                [0.0, -0.86602500000000004, 0.5],
                [-0.22414400000000001, -0.83651600000000004, 0.5],
                [-0.43301299999999998, -0.75, 0.5],
                [-0.61237200000000003, -0.61237200000000003, 0.5],
                [-0.75, -0.43301299999999998, 0.5],
                [-0.83651600000000004, -0.22414400000000001, 0.5],
                [-0.86602500000000004, 0.0, 0.5],
                [-0.83651600000000004, 0.22414400000000001, 0.5],
                [-0.75, 0.43301299999999998, 0.5],
                [-0.61237200000000003, 0.61237200000000003, 0.5],
                [-0.43301299999999998, 0.75, 0.5],
                [-0.22414400000000001, 0.83651600000000004, 0.5],
                [0.0, 0.74314499999999994, 0.66913100000000003],
                [0.19234000000000001, 0.71782299999999999, 0.66913100000000003],
                [0.37157200000000001, 0.64358199999999999, 0.66913100000000003],
                [0.52548300000000003, 0.52548300000000003, 0.66913100000000003],
                [0.64358199999999999, 0.37157200000000001, 0.66913100000000003],
                [0.71782299999999999, 0.19234000000000001, 0.66913100000000003],
                [0.74314499999999994, 0.0, 0.66913100000000003],
                [0.71782299999999999, -0.19234000000000001, 0.66913100000000003],
                [0.64358199999999999, -0.37157200000000001, 0.66913100000000003],
                [0.52548300000000003, -0.52548300000000003, 0.66913100000000003],
                [0.37157200000000001, -0.64358199999999999, 0.66913100000000003],
                [0.19234000000000001, -0.71782299999999999, 0.66913100000000003],
                [0.0, -0.74314499999999994, 0.66913100000000003],
                [-0.19234000000000001, -0.71782299999999999, 0.66913100000000003],
                [-0.37157200000000001, -0.64358199999999999, 0.66913100000000003],
                [-0.52548300000000003, -0.52548300000000003, 0.66913100000000003],
                [-0.64358199999999999, -0.37157200000000001, 0.66913100000000003],
                [-0.71782299999999999, -0.19234000000000001, 0.66913100000000003],
                [-0.74314499999999994, 0.0, 0.66913100000000003],
                [-0.71782299999999999, 0.19234000000000001, 0.66913100000000003],
                [-0.64358199999999999, 0.37157200000000001, 0.66913100000000003],
                [-0.52548300000000003, 0.52548300000000003, 0.66913100000000003],
                [-0.37157200000000001, 0.64358199999999999, 0.66913100000000003],
                [-0.19234000000000001, 0.71782299999999999, 0.66913100000000003],
                [0.0, 0.587785, 0.80901699999999999],
                [0.20103399999999999, 0.55233699999999997, 0.80901699999999999],
                [0.37782100000000002, 0.45027, 0.80901699999999999],
                [0.50903699999999996, 0.29389300000000002, 0.80901699999999999],
                [0.57885500000000001, 0.10206800000000001, 0.80901699999999999],
                [0.57885500000000001, -0.10206800000000001, 0.80901699999999999],
                [0.50903699999999996, -0.29389300000000002, 0.80901699999999999],
                [0.37782100000000002, -0.45027, 0.80901699999999999],
                [0.20103399999999999, -0.55233699999999997, 0.80901699999999999],
                [0.0, -0.587785, 0.80901699999999999],
                [-0.20103399999999999, -0.55233699999999997, 0.80901699999999999],
                [-0.37782100000000002, -0.45027, 0.80901699999999999],
                [-0.50903699999999996, -0.29389300000000002, 0.80901699999999999],
                [-0.57885500000000001, -0.10206800000000001, 0.80901699999999999],
                [-0.57885500000000001, 0.10206800000000001, 0.80901699999999999],
                [-0.50903699999999996, 0.29389300000000002, 0.80901699999999999],
                [-0.37782100000000002, 0.45027, 0.80901699999999999],
                [-0.20103399999999999, 0.55233699999999997, 0.80901699999999999],
                [0.0, 0.40673700000000002, 0.91354500000000005],
                [0.20336799999999999, 0.352244, 0.91354500000000005],
                [0.352244, 0.20336799999999999, 0.91354500000000005],
                [0.40673700000000002, 0.0, 0.91354500000000005],
                [0.352244, -0.20336799999999999, 0.91354500000000005],
                [0.20336799999999999, -0.352244, 0.91354500000000005],
                [0.0, -0.40673700000000002, 0.91354500000000005],
                [-0.20336799999999999, -0.352244, 0.91354500000000005],
                [-0.352244, -0.20336799999999999, 0.91354500000000005],
                [-0.40673700000000002, 0.0, 0.91354500000000005],
                [-0.352244, 0.20336799999999999, 0.91354500000000005],
                [-0.20336799999999999, 0.352244, 0.91354500000000005],
                [0.0, 0.20791200000000001, 0.97814800000000002],
                [0.18005699999999999, 0.10395600000000001, 0.97814800000000002],
                [0.18005699999999999, -0.10395600000000001, 0.97814800000000002],
                [0.0, -0.20791200000000001, 0.97814800000000002],
                [-0.18005699999999999, -0.10395600000000001, 0.97814800000000002],
                [-0.18005699999999999, 0.10395600000000001, 0.97814800000000002],
                [0.0, 0.0, 1]
            ])
            self.patch_count = 145
            self.row_conversion_factor = np.array([
                0.0435449227,
                0.0416418006,
                0.0473984151,
                0.0406730411,
                0.0428934136,
                0.0445221864,
                0.0455168385,
                0.0344199465
            ])
            self.row_patches = np.array([
                30,
                30,
                24,
                24,
                18,
                12,
                6,
                1,
            ])
            self.patch_conversion_factor = np.array([
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.04354492, 0.04354492, 0.04354492, 0.04354492, 0.04354492,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.0416418, 0.0416418, 0.0416418, 0.0416418, 0.0416418,
                0.04739842, 0.04739842, 0.04739842, 0.04739842, 0.04739842,
                0.04739842, 0.04739842, 0.04739842, 0.04739842, 0.04739842,
                0.04739842, 0.04739842, 0.04739842, 0.04739842, 0.04739842,
                0.04739842, 0.04739842, 0.04739842, 0.04739842, 0.04739842,
                0.04739842, 0.04739842, 0.04739842, 0.04739842, 0.04067304,
                0.04067304, 0.04067304, 0.04067304, 0.04067304, 0.04067304,
                0.04067304, 0.04067304, 0.04067304, 0.04067304, 0.04067304,
                0.04067304, 0.04067304, 0.04067304, 0.04067304, 0.04067304,
                0.04067304, 0.04067304, 0.04067304, 0.04067304, 0.04067304,
                0.04067304, 0.04067304, 0.04067304, 0.04289341, 0.04289341,
                0.04289341, 0.04289341, 0.04289341, 0.04289341, 0.04289341,
                0.04289341, 0.04289341, 0.04289341, 0.04289341, 0.04289341,
                0.04289341, 0.04289341, 0.04289341, 0.04289341, 0.04289341,
                0.04289341, 0.04452219, 0.04452219, 0.04452219, 0.04452219,
                0.04452219, 0.04452219, 0.04452219, 0.04452219, 0.04452219,
                0.04452219, 0.04452219, 0.04452219, 0.04551684, 0.04551684,
                0.04551684, 0.04551684, 0.04551684, 0.04551684, 0.03441995
            ])
        # Set radiation rose metrics
        self.radiation_rose_vectors = np.array([
            [0.0, 1.0, 0.0],
            [-0.17364817766693033, 0.98480775301220802, 0.0],
            [-0.34202014332566871, 0.93969262078590843, 0.0],
            [-0.49999999999999994, 0.86602540378443871, 0.0],
            [-0.64278760968653925, 0.76604444311897801, 0.0],
            [-0.76604444311897801, 0.64278760968653936, 0.0],
            [-0.8660254037844386, 0.50000000000000011, 0.0],
            [-0.93969262078590832, 0.34202014332566882, 0.0],
            [-0.98480775301220802, 0.17364817766693041, 0.0],
            [-1.0, 0.0, 0.0],
            [-0.98480775301220802, -0.1736481776669303, 0.0],
            [-0.93969262078590843, -0.34202014332566871, 0.0],
            [-0.86602540378443871, -0.49999999999999978, 0.0],
            [-0.76604444311897801, -0.64278760968653936, 0.0],
            [-0.64278760968653947, -0.7660444431189779, 0.0],
            [-0.49999999999999994, -0.86602540378443871, 0.0],
            [-0.34202014332566888, -0.93969262078590832, 0.0],
            [-0.17364817766693028, -0.98480775301220802, 0.0],
            [0.0, -1.0, 0.0],
            [0.17364817766693047, -0.98480775301220802, 0.0],
            [0.34202014332566866, -0.93969262078590843, 0.0],
            [0.50000000000000011, -0.8660254037844386, 0.0],
            [0.64278760968653925, -0.76604444311897801, 0.0],
            [0.7660444431189779, -0.64278760968653947, 0.0],
            [0.86602540378443837, -0.50000000000000044, 0.0],
            [0.93969262078590821, -0.34202014332566938, 0.0],
            [0.98480775301220802, -0.17364817766693033, 0.0],
            [1.0, 0.0, 0.0],
            [0.98480775301220813, 0.17364817766692997, 0.0],
            [0.93969262078590854, 0.34202014332566816, 0.0],
            [0.8660254037844386, 0.50000000000000011, 0.0],
            [0.76604444311897812, 0.64278760968653925, 0.0],
            [0.64278760968653958, 0.76604444311897779, 0.0],
            [0.50000000000000044, 0.86602540378443837, 0.0],
            [0.3420201433256686, 0.93969262078590843, 0.0],
            [0.17364817766693127, 0.98480775301220791, 0.0],
        ])
        self.radiation_rose_angles = np.radians(np.arange(0, 360, 10))[::-1]
        # Migrate the sky dome stuff up here - methods to add next
        self.direct_sky_matrix = self.gendaymtx(direct=True)
        print("Direct sky matrix calculated: {0:}".format(pathlib.Path(self.wea_file).with_suffix(".dirmtx")))
        self.diffuse_sky_matrix = self.gendaymtx(direct=False)
        print("Diffuse sky matrix calculated: {0:}".format(pathlib.Path(self.wea_file).with_suffix(".diffmtx")))
        self.total_sky_matrix = self.direct_sky_matrix + self.diffuse_sky_matrix

        return self

    # Methods below here for NV method

    

    # Methods below here are for plotting only

    def plot_heatmap(self, variable, cmap='Greys', close=True, save=False):
        # Construct the save_path and create directory if it doesn't exist
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "heatmap_{}.png".format(variable)

        # Instantiate figure
        fig, ax = plt.subplots(1, 1, figsize=(15, 5))

        # Load data
        series = self.df[variable].to_frame()

        # Remove timezone-awareness
        try:
            series.index = series.index.tz_convert(None)
        except Exception as e:
            raise e

        # Reshape data into time/day matrix
        ll = series.pivot_table(columns=series.index.date, index=series.index.time).values[::-1]

        # Plot data
        heatmap = ax.imshow(
            ll,
            extent=[mdates.date2num(series.index.min()), mdates.date2num(series.index.max()), 726449, 726450],
            aspect='auto',
            cmap=cmap,
            interpolation='none'
        )

        # Formatting
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.yaxis_date()
        ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.invert_yaxis()
        ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
        plt.setp(ax.get_xticklabels(), ha='left', color='k')
        plt.setp(ax.get_yticklabels(), color='k')
        [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]
        ax.grid(b=True, which='major', color='white', linestyle='-', alpha=1)
        cb = fig.colorbar(heatmap, orientation='horizontal', drawedges=False, fraction=0.05, aspect=100, pad=0.075)
        plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='k')
        cb.outline.set_visible(False)
        plt.title("{}\n{} - {} - {}".format(renamer(series.columns[0]), self.city, self.country, self.station_id), color='k', y=1.01)
        plt.tight_layout()

        # Save figure
        if save:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
            print("Heatmap saved to {}".format(save_path))
        if close:
            plt.close()

        return fig

    def plot_diurnal(self, dew_point=False, close=True, save=False):
        # Construct the save_path and create directory if it doesn't exist
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "diurnal_{}.png".format("dpt" if dew_point else "rh")

        # Group dry-bulb temperatures
        a_gp = self.dry_bulb_temperature.groupby([self.index.month, self.index.hour])
        a_min = a_gp.min().reset_index(drop=True)
        a_mean = a_gp.mean().reset_index(drop=True)
        a_max = a_gp.max().reset_index(drop=True)

        # Group relative humidity / dewpoint temperature
        if dew_point:
            b_gp = self.dew_point_temperature.groupby([self.index.month, self.index.hour])
            b_var = "Dew-point Temperature (°C)"
        else:
            b_gp = self.relative_humidity.groupby([self.relative_humidity.index.month, self.relative_humidity.index.hour])
            b_var = 'Relative Humidity (%)'
        b_min = b_gp.min().reset_index(drop=True)
        b_mean = b_gp.mean().reset_index(drop=True)
        b_max = b_gp.max().reset_index(drop=True)

        # Group solar radiation
        c_global_mean = self.global_horizontal_radiation.groupby([self.index.month, self.index.hour]).mean().reset_index(drop=True)
        c_diffuse_mean = self.diffuse_horizontal_radiation.groupby([self.index.month, self.index.hour]).mean().reset_index(drop=True)
        c_direct_mean = self.direct_normal_radiation.groupby([self.index.month, self.index.hour]).mean().reset_index(drop=True)

        # Instantiate plot
        fig, ax = plt.subplots(3, 1, figsize=(15, 8))

        # Plot DBT
        [ax[0].plot(a_mean.iloc[i:i + 24], color='#BC204B', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
        [ax[0].fill_between(np.arange(i, i + 24), a_min.iloc[i:i + 24], a_max.iloc[i:i + 24], color='#BC204B', alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
        ax[0].set_ylabel('Dry-bulb Temperature (°C)', labelpad=2, color='k')
        ax[0].yaxis.set_major_locator(MaxNLocator(7))

        # Plot DPT / RH
        [ax[1].plot(b_mean.iloc[i:i + 24], color='#00617F', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
        [ax[1].fill_between(np.arange(i, i + 24), b_min.iloc[i:i + 24], b_max.iloc[i:i + 24], color='#00617F', alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
        ax[1].set_ylabel(b_var, labelpad=2, color='k')
        ax[1].yaxis.set_major_locator(MaxNLocator(7))
        if not dew_point:
            ax[1].set_ylim([0, 100])

        # Plot solar
        [ax[2].plot(c_direct_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls='--', label='Direct Normal Radiation') for i in
         np.arange(0, 288)[::24]]
        [ax[2].plot(c_diffuse_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls=':', label='Diffuse Horizontal Radiation') for i
         in np.arange(0, 288)[::24]]
        [ax[2].plot(c_global_mean.iloc[i:i + 24], color='#FF8F1C', lw=2, ls='-', label='Global Horizontal Radiation') for
         i in np.arange(0, 288)[::24]]
        ax[2].set_ylabel('Solar Radiation (W/m$^{2}$)', labelpad=2, color='k')
        ax[2].yaxis.set_major_locator(MaxNLocator(7))

        # Format plot area
        [[i.spines[spine].set_visible(False) for spine in ['top', 'right']] for i in ax]
        [[i.spines[j].set_color('k') for i in ax] for j in ['bottom', 'left']]
        [i.xaxis.set_ticks(np.arange(0, 288, 24)) for i in ax]
        [i.set_xlim([0, 287]) for i in ax]
        [plt.setp(i.get_yticklabels(), color='k') for i in ax]
        [i.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
                           ha='left', color='k') for i in ax]
        [i.get_xaxis().set_ticklabels([]) for i in [ax[0], ax[1]]]
        [i.grid(b=True, which='major', axis='both', c='k', ls='--', lw=1, alpha=0.3) for i in ax]
        [i.tick_params(length=0) for i in ax]

        ax[2].set_ylim([0, ax[2].get_ylim()[1]])

        # Legend
        handles, labels = ax[2].get_legend_handles_labels()
        lgd = ax[2].legend(bbox_to_anchor=(0.5, -0.2), loc=8, ncol=3, borderaxespad=0., frameon=False, handles=[handles[0], handles[12], handles[24]], labels=[labels[0], labels[12], labels[24]])
        lgd.get_frame().set_facecolor((1, 1, 1, 0))
        [plt.setp(text, color='k') for text in lgd.get_texts()]

        # Add a title
        title = plt.suptitle(
            "Monthly average diurnal profile\n{0:} - {1:} - {2:}".format(self.city, self.country, self.station_id),
            color='k', y=1.025)

        # Tidy plot
        plt.tight_layout()

        # Save figure
        if save:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
            print("Diurnal plot saved to {}".format(save_path))
        if close:
            plt.close()

        return fig

    def plot_radiation_rose(self, save=False, close=False):

        direct_values = self.direct_sky_matrix.sum(axis=0)
        diffuse_values = self.diffuse_sky_matrix.sum(axis=0)

        direct_sky_radiation_rose_values = []
        diffuse_sky_radiation_rose_values = []
        for vec in self.radiation_rose_vectors:
            direct_radiation = 0
            diffuse_radiation = 0
            for patch_number, patch_vector in enumerate(self.patch_vectors):
                vector_angle = angle_between(patch_vector, vec, degrees=True)
                if vector_angle < 90:
                    direct_radiation += direct_values[patch_number] * np.cos(np.radians(vector_angle))
                    diffuse_radiation += diffuse_values[patch_number] * np.cos(np.radians(vector_angle))
            direct_sky_radiation_rose_values.append(direct_radiation)
            diffuse_sky_radiation_rose_values.append(diffuse_radiation)
        direct_sky_radiation_rose_values = np.array(direct_sky_radiation_rose_values)
        diffuse_sky_radiation_rose_values = np.array(diffuse_sky_radiation_rose_values)
        total_sky_radiation_rose_values = direct_sky_radiation_rose_values + diffuse_sky_radiation_rose_values

        max_val = total_sky_radiation_rose_values.max() / 1000
        figs = []
        for tx, vals in {"Direct Radiation": direct_sky_radiation_rose_values, "Diffuse Radiation": diffuse_sky_radiation_rose_values, "Total Radiation": total_sky_radiation_rose_values}.items():

            # Construct the save_path and create directory if it doesn't exist
            save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "radiationrose_{}.png".format(tx)

            vals = vals / 1000
            colors = [cm.OrRd(i) for i in np.interp(vals, [min(vals), max(vals)], [0, 1])]
            fig, ax = plt.subplots(1, 1, figsize=(6, 6), subplot_kw={'projection': "polar"})
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            ax.bar(
                self.radiation_rose_angles,
                vals,
                width=(np.pi * 2 / 36), zorder=5, bottom=0.0, color=colors, alpha=1, edgecolor="w", linewidth=0)
            ax.set_ylim(0, max_val)
            ax.spines['polar'].set_visible(False)
            ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
            ti = ax.set_title("{}\n{} - {} - {}".format(tx, self.city, self.country, self.station_id), color="k", loc="left", va="bottom", ha="left", fontsize="large", y=1)
            plt.tight_layout()

            # Save figure
            if save:
                save_path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
                print("Radiation rose saved to {}".format(save_path))
            if close:
                plt.close()
            figs.append(fig)

        return figs

    def plot_windrose(self, season_period="Annual", day_period="Daily", nsector=16, cmap=None, save=False, close=False):
        from windrose import WindroseAxes

        # Construct the save_path and create directory if it doesn't exist
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "windrose_{}_{}.png".format(season_period, day_period)
        
        # Descibe a set of masks to remove unwanted hours of the year
        toy_masks = {
            "Daily": ((self.index.hour >= 0) & (self.index.hour <= 24)),
            "Morning": ((self.index.hour >= 5) & (self.index.hour <= 10)),
            "Midday": ((self.index.hour >= 11) & (self.index.hour <= 13)),
            "Afternoon": ((self.index.hour >= 14) & (self.index.hour <= 18)),
            "Evening": ((self.index.hour >= 19) & (self.index.hour <= 22)),
            "Night": ((self.index.hour >= 23) | (self.index.hour <= 4)),
    
            "Annual": ((self.index.month >= 1) & (self.index.month <= 12)),
            "Spring": ((self.index.month >= 3) & (self.index.month <= 5)),
            "Summer": ((self.index.month >= 6) & (self.index.month <= 8)),
            "Autumn": ((self.index.month >= 9) & (self.index.month <= 11)),
            "Winter": ((self.index.month <= 2) | (self.index.month >= 12))
        }
        speed_mask = (self.wind_speed != 0)
        direction_mask = (self.wind_direction != 0)
        mask = np.array([toy_masks[day_period], toy_masks[season_period], speed_mask, direction_mask]).all(axis=0)
    
        fig = plt.figure(figsize=(6, 6))
        ax = WindroseAxes.from_ax()
        unit = "m/s"
        ax.bar(self.wind_direction[mask], self.wind_speed[mask], normed=True,
               bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], opening=1, edgecolor='White',
               lw=0.25, nsector=nsector,
               cmap=plt.cm.Purples if cmap is None else cmap)
    
        lgd = ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', frameon=False, title=unit)
        lgd.get_frame().set_facecolor((1, 1, 1, 0))
        [plt.setp(text, color='k') for text in lgd.get_texts()]
        plt.setp(lgd.get_title(), color='k')
    
        for i, leg in enumerate(lgd.get_texts()):
            b = leg.get_text().replace('[', '').replace(')', '').split(' : ')
            lgd.get_texts()[i].set_text(b[0] + ' to ' + b[1])
    
        ax.grid(linestyle=':', color='k', alpha=0.5)
        ax.spines['polar'].set_visible(False)
        plt.setp(ax.get_xticklabels(), color='k')
        plt.setp(ax.get_yticklabels(), color='k')
        ax.set_title("{2:} - {3:} - {4:}\n{0:} - {1:} - {5:}".format(self.city, self.country, season_period, day_period, "Wind speed", self.station_id), y=1.06, color='k')
    
        plt.tight_layout()

        # Save figure
        if save:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
            print("Windrose saved to {}".format(save_path))
        if close:
            plt.close()

        return fig

    # TODO: plot_wind_weibull, plot_utci_frequency, plot_utci_heatmap
