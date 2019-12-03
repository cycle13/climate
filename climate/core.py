import os
import pathlib
import subprocess
from io import StringIO
import pandas as pd
import numpy as np

from .helpers import *


class BuroHappold(object):

    def __init__(self):

        self.colors = {
            "official": {
                "BrightGreen": [108, 194, 78],
                "BrightPurple": [112, 47, 138],
                "BrightRed": [213, 0, 50],
                "DarkBlue": [0, 97, 127],
                "DarkRed": [126, 45, 64],
                "Green": [34, 136, 72],
                "Grey": [150, 140, 131],
                "LightBlue": [141, 185, 202],
                "LightRed": [188, 32, 75],
                "Orange": [255, 143, 28],
                "Pink": [231, 130, 169],
                "Purple": [36, 19, 95],
                "BrandGreen": [196, 214, 0],
                "Teal": [0, 164, 153],
            },
            "ppt": {
                "Purple": [143, 114, 176],
                "BrandGreen": [196, 214, 0],
                "Pink": [230, 49, 135],
                "LightBlue": [0, 176, 240],
                "DarkBlue": [28, 54, 96],
                "Orange": [235, 103, 28],
                "DarkPurple": [85, 62, 112],
                "DarkGreen": [118, 128, 0],
                "DarkPink": [149, 18, 80],
                "Blue": [0, 106, 144],
                "ReallyDarkBlue": [17, 32, 58],
                "Brown": [145, 61, 13],
            },
        }


class EPW(object):

    def __init__(self, file_path):

        # Metadata
        self.filepath = os.path.abspath(file_path) if file_path is not None else None

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

        # Derived variables
        self.solar_altitude = None
        self.solar_azimuth = None
        self.wet_bulb_temperature = None
        self.enthalpy = None
        self.humidity_ratio = None
        self.universal_thermal_climate_index = None
        self.standard_effective_temperature = None
        self.solar_adjusted_mean_radiant_temperature = None
        self.mean_radiant_temperature = None
        self.ground_temperature_at_depth = None
        self.direct_sky_matrix = None
        self.diffuse_sky_matrix = None
        self.total_sky_matrix = None

        # Path variables
        self.wea_file = None
        self.csv_file = None

    def read(self):
        with open(self.filepath, "r") as f:
            dat = f.readlines()

            # Read location data
            self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = dat[0].strip().split(",")[1:]
            self.latitude = float(self.latitude)
            self.longitude = float(self.longitude)
            self.time_zone = float(self.time_zone)
            self.elevation = float(self.elevation)

            # Read design conditions data
            self.design_conditions = ",".join(dat[1].strip().split(",")[1:])

            # Read typical extreme periods data
            self.typical_extreme_periods = ",".join(dat[2].strip().split(",")[1:])

            # Read ground temperatures data from weatherfile
            self.ground_temperatures = ",".join(dat[3].strip().split(",")[1:])

            # Read holidays/daylight savings data
            self.holidays_daylight_savings = ",".join(dat[4].strip().split(",")[1:])

            # Read comments 1 data
            self.comments_1 = ",".join(dat[5].strip().split(",")[1:])

            # Read comments 2 data
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
            df.index = pd.date_range(
                start="2018-01-01 00:00:00",
                end="2019-01-01 00:00:00",
                freq="60T",
                closed="left").tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60))

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

    def to_csv(self, file_path=None):
        if file_path is None:
            file_path = pathlib.Path(self.filepath).with_suffix(".csv")
        self.df.to_csv(file_path)
        self.csv_file = str(file_path)
        print("CSV file created: {0:}".format(self.csv_file))
        return self.csv_file

    def to_wea(self, file_path=None):
        if file_path is None:
            file_path = pathlib.Path(self.filepath).with_suffix(".wea")
        header = "place {0:}_{1:}\nlatitude {2:0.2f}\nlongitude {3:0.2f}\ntime_zone {4:0.0f}\nsite_elevation {5:0.2f}\nweather_data_file_units 1".format(slugify(self.city), slugify(self.country), self.latitude, self.longitude, -self.time_zone / 15, self.elevation)
        values = []
        for n, dt in enumerate(self.df.index):
            values.append("{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n], self.diffuse_horizontal_radiation[n]))
        with open(file_path, "w") as f:
            f.write(header + "\n" + "\n".join(values) + "\n")
        self.wea_file = str(file_path)
        print("WEA file created: {0:}".format(self.wea_file))
        return self.wea_file

    def run_psychrometrics(self):
        self.wet_bulb_temperature = self.df.apply(lambda x: humidity_ratio_relative_humidity(dry_bulb_temperature=x["dry_bulb_temperature"], relative_humidity=x["relative_humidity"] / 100, pressure=x["atmospheric_station_pressure"] / 1000), axis=1)
        print("Wet-bulb temperature calculated")
        self.enthalpy = self.df.apply(lambda x: enthalpy_relative_humidity(dry_bulb_temperature=x["dry_bulb_temperature"], relative_humidity=x["relative_humidity"] / 100, pressure=x["atmospheric_station_pressure"] / 1000), axis=1)
        print("Enthalpy calculated")
        self.humidity_ratio = self.df.apply(lambda x: wet_bulb_temperature_relative_humidity(dry_bulb_temperature=x["dry_bulb_temperature"], relative_humidity=x["relative_humidity"] / 100, pressure=x["atmospheric_station_pressure"] / 1000), axis=1)
        print("Humidity ratio calculated")
        self.wet_bulb_temperature.name = "wet_bulb_temperature"
        self.enthalpy.name = "enthalpy"
        self.humidity_ratio.name = "humidity_ratio"
        self.df = pd.concat([self.df, self.wet_bulb_temperature, self.enthalpy, self.humidity_ratio], axis=1)

    def run_sunposition(self):
        from pysolar.solar import get_altitude_fast, get_azimuth_fast
        self.solar_altitude = pd.Series(index=self.df.index, data=[float(get_altitude_fast(self.latitude, self.longitude, i)) for i in self.df.index], name="solar_altitude")
        print("Solar altitude calculated")
        self.solar_azimuth = pd.Series(index=self.df.index, data=[float(get_azimuth_fast(self.latitude, self.longitude, i)) for i in self.df.index], name="solar_azimuth")
        print("Solar azimuth calculated")
        self.df = pd.concat([self.df, self.solar_altitude, self.solar_azimuth], axis=1)

    def run_ground_temperatures(self, depth=0.1):
        mn = self.dry_bulb_temperature.mean()
        rng = self.dry_bulb_temperature.max() - self.dry_bulb_temperature.min()
        coldest_day = self.dry_bulb_temperature.resample("1D").mean().idxmin()

        gndts = []
        for j in [i if i > 0 else i + 365 for i in (self.df.index - coldest_day).total_seconds() / 86400]:
            gndts.append(ground_temperature_at_depth(depth, mn, rng, j, soil_diffusivity=0.01))
        self.ground_temperature_at_depth = pd.Series(name="ground_temperature_at_depth", index=self.df.index, data=gndts)
        print("Ground temperature calculated")
        self.df = pd.concat([self.df, self.ground_temperature_at_depth], axis=1)

    def run_pedestrian_wind_speed(self):
        self.pedestrian_wind_speed = pd.Series(name="pedestrian_wind_speed", index=self.df.index, data=[wind_speed_at_height(ws=i, h1=10, h2=1.5) for i in self.wind_speed])
        print("Pedestrian wind-speed calculated")
        self.df = pd.concat([self.df, self.pedestrian_wind_speed], axis=1)

    def run_gendaymtx(self, reinhart=True):

        if self.wea_file is None:
            self.to_wea()

        diff_mtx_file = pathlib.Path(self.wea_file).with_suffix(".diffmtx")
        dir_mtx_file = pathlib.Path(self.wea_file).with_suffix(".dirmtx")

        # Create command strings
        diffuse_mtx_cmd = '"C:/Radiance/bin/gendaymtx" -m {2:} -s -O1 "{0:}" > "{1:}"'.format(self.wea_file, diff_mtx_file,
                                                                                              2 if reinhart else 1)
        direct_mtx_cmd = '"C:/Radiance/bin/gendaymtx" -m {2:} -d -O1 "{0:}" > "{1:}"'.format(self.wea_file, dir_mtx_file,
                                                                                             2 if reinhart else 1)

        # Run commands
        subprocess.call(diffuse_mtx_cmd, shell=True)
        subprocess.call(direct_mtx_cmd, shell=True)

        # Load matrices
        num_of_patches_in_each_row = {
            1: np.array([30, 30, 24, 24, 18, 12, 6, 1]),
            2: np.array([60, 60, 60, 60, 48, 48, 48, 48, 36, 36, 24, 24, 12, 12, 1])
        }
        patch_conversion_factor = {
            1: np.array(
                [0.0435449227, 0.0416418006, 0.0473984151, 0.0406730411, 0.0428934136, 0.0445221864, 0.0455168385,
                 0.0344199465]),
            2: np.array(
                [0.0113221971, 0.0111894547, 0.0109255262, 0.0105335058, 0.0125224872, 0.0117312774, 0.0108025291,
                 0.00974713106, 0.011436609, 0.00974295956, 0.0119026242, 0.00905126163, 0.0121875626, 0.00612971396,
                 0.00921483254])
        }

        # Create the conversion factor for each patch
        if reinhart:
            conversion_factor = np.repeat(patch_conversion_factor[2], num_of_patches_in_each_row[2])
        else:
            conversion_factor = np.repeat(patch_conversion_factor[1], num_of_patches_in_each_row[1])

        with open(diff_mtx_file, "r") as f:
            d = [i.strip() for i in f.readlines()][8:]
        diff_chunks = np.array([[j.split(" ") for j in i[:-1]] for i in chunk(d, 8761)])[1:].astype(np.float64)
        diff_chunks = np.sum(np.multiply(diff_chunks, [0.265074126, 0.670114631, 0.064811243]), axis=2)
        diff_chunks = np.multiply(diff_chunks.T, conversion_factor)
        print("Diffuse radiation matrix created: {0:}".format(str(diff_mtx_file)))

        with open(dir_mtx_file, "r") as f:
            d = [i.strip() for i in f.readlines()][8:]
        dir_chunks = np.array([[j.split(" ") for j in i[:-1]] for i in chunk(d, 8761)])[1:].astype(np.float64)
        dir_chunks = np.sum(np.multiply(dir_chunks, [0.265074126, 0.670114631, 0.064811243]), axis=2)
        dir_chunks = np.multiply(dir_chunks.T, conversion_factor)
        print("Direct radiation matrix created: {0:}".format(str(dir_mtx_file)))

        self.direct_sky_matrix = dir_chunks
        self.diffuse_sky_matrix = diff_chunks
        self.total_sky_matrix = dir_chunks + diff_chunks

    # def
