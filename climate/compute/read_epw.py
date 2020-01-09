import pandas as pd
from io import StringIO

from climate.common.constants import DATETIME_INDEX

def read_epw_normal(self):
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
    with open(self.file_path.absolute(), "r") as f:
        data = f.readlines()

        # Read location data
        self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = \
        data[0].strip().split(",")[1:]
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)
        self.time_zone = float(self.time_zone)
        self.elevation = float(self.elevation)
        self.design_conditions = ",".join(data[1].strip().split(",")[1:])
        self.typical_extreme_periods = ",".join(data[2].strip().split(",")[1:])
        self.ground_temperatures = ",".join(data[3].strip().split(",")[1:])
        self.holidays_daylight_savings = ",".join(data[4].strip().split(",")[1:])
        self.comments_1 = ",".join(data[5].strip().split(",")[1:])
        self.comments_2 = ",".join(data[6].strip().split(",")[1:])

        # Read the data table
        df = pd.read_csv(StringIO("\n".join(data[8:])), header=None)

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
        self.index = DATETIME_INDEX.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(
            hours=self.time_zone)
        df.index = self.index

        # Drop date/time columns
        df.drop(columns=["year", "month", "day", "hour", "minute"], inplace=True)

        # Make loaded data accessible
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
    return self

def read_epw_ccwwg(self):
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
        data = f.readlines()

        # Read location data
        self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = \
        data[0].strip().split(",")[1:]
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)
        self.time_zone = float(self.time_zone)
        self.elevation = float(self.elevation)
        self.design_conditions = ",".join(data[1].strip().split(",")[1:])
        self.typical_extreme_periods = ",".join(data[2].strip().split(",")[1:])
        self.ground_temperatures = ",".join(data[3].strip().split(",")[1:])
        self.holidays_daylight_savings = ",".join(data[4].strip().split(",")[1:])
        self.comments_1 = ",".join(data[5].strip().split(",")[1:])
        self.comments_2 = ",".join(data[6].strip().split(",")[1:])

        # Read the data table
        df = pd.read_csv(StringIO("\n".join(data[8:])), header=None)

        print(df.shape)

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
        ]

        # Create datetime index - using 2018 as base year (a Monday starting year without leap-day)
        self.index = DATETIME_INDEX.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(
            hours=self.time_zone)
        df.index = self.index

        # Drop date/time columns
        df.drop(columns=["year", "month", "day", "hour", "minute"], inplace=True)

        # Make loaded data accessible
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
    return self

def read_epw_cibse(file_path):
    return None
