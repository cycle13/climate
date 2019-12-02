import os
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

    def __init__(self, filepath):

        # Metadata
        self.filepath = os.path.abspath(filepath) if filepath is not None else None

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
        self.solar_wet_bulb_temperature = None
        self.solar_enthalpy = None
        self.solar_humidity_ratio = None
        self.universal_thermal_climate_index = None
        self.standard_effective_temperature = None
        self.solar_adjusted_mean_radiant_temperature = None
        self.mean_radiant_temperature = None


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

            # Read ground temperatures data
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


    def to_csv(self, filepath):
        self.df.to_csv(filepath)
        return filepath
