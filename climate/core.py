import os
from io import StringIO
import subprocess

import pandas as pd
import numpy as np

from .helpers import *


class EPW(object):

    def __init__(self, filepath):

        # Metadata
        self.filepath = os.path.normpath(filepath) if filepath is not None else None

        # Header information
        self.city = ''
        self.region = ''
        self.country = ''
        self.dataset_type = ''
        self.station_id = ''
        self.latitude = []
        self.longitude = []
        self.elevation = 0
        self.time_zone = 0
        self.design_conditions = ''
        self.typical_extreme_periods = ''
        self.ground_temperatures = ''
        self.holidays_daylight_savings = ''
        self.comments_1 = ''
        self.comments_2 = ''

        # Variables information
        self.year = []
        self.month = []
        self.day = []
        self.hour = []
        self.minute = []
        self.data = None
        self.data_source_and_uncertainty_flags = []
        self.dry_bulb_temperature = []
        self.dew_point_temperature = []
        self.relative_humidity = []
        self.atmospheric_station_pressure = []
        self.extraterrestrial_horizontal_radiation = []
        self.extraterrestrial_direct_normal_radiation = []
        self.horizontal_infrared_radiation_intensity = []
        self.global_horizontal_radiation = []
        self.direct_normal_radiation = []
        self.diffuse_horizontal_radiation = []
        self.global_horizontal_illuminance = []
        self.direct_normal_illuminance = []
        self.diffuse_horizontal_illuminance = []
        self.zenith_luminance = []
        self.wind_direction = []
        self.wind_speed = []
        self.total_sky_cover = []
        self.opaque_sky_cover = []
        self.visibility = []
        self.ceiling_height = []
        self.present_weather_observation = []
        self.present_weather_codes = []
        self.precipitable_water = []
        self.aerosol_optical_depth = []
        self.snow_depth = []
        self.days_since_last_snowfall = []
        self.albedo = []
        self.liquid_precipitation_depth = []
        self.liquid_precipitation_quantity = []

    def read(self):
        with open(self.filepath, "r") as f:
            dat = f.readlines()

            # Read location data
            self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = dat[0].strip().split(",")[1:]
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
            df.columns = ['year', 'month', 'day', 'hour', 'minute', 'data_source_and_uncertainty_flags',
                          'dry_bulb_temperature', 'dew_point_temperature', 'relative_humidity',
                          'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
                          'extraterrestrial_direct_normal_radiation', 'horizontal_infrared_radiation_intensity',
                          'global_horizontal_radiation', 'direct_normal_radiation', 'diffuse_horizontal_radiation',
                          'global_horizontal_illuminance', 'direct_normal_illuminance',
                          'diffuse_horizontal_illuminance', 'zenith_luminance', 'wind_direction', 'wind_speed',
                          'total_sky_cover', 'opaque_sky_cover', 'visibility', 'ceiling_height',
                          'present_weather_observation', 'present_weather_codes', 'precipitable_water',
                          'aerosol_optical_depth', 'snow_depth', 'days_since_last_snowfall', 'albedo',
                          'liquid_precipitation_depth', 'liquid_precipitation_quantity', ]
            # Create datetime index - using 2018 as base year (a Monday starting year without leap-day)
            df.index = pd.date_range(start="2018-01-01 00:00:00", end="2019-01-01 00:00:00", freq="60T", closed="left",
                                     tz=int(self.time_zone * 60 * 60))
            self.dt = df.index
            df.drop(columns=["year", "month", "day", "hour", "minute"], inplace=True)

            # Make loaded data accesible
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

