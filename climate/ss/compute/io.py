import pandas as pd
from io import StringIO
import pathlib
import pickle
import warnings

# from climate.core import Weather
from climate import slugify
from climate import DATETIME_INDEX
from climate import weatherfile_ground_temperatures
from climate import wind_speed_at_height


def load_epw(self, kind="standard"):
    """
        Read EPW weather-file into weather object.

        Parameters
        ----------


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

        # Create datetime index - using 2018 as base year (a Monday starting year without leap-day)
        self.index = DATETIME_INDEX.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(
            hours=self.time_zone)
        df.index = self.index

        if kind == "standard":

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

        elif kind == "ccwwg":

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

        elif kind == "cibse":

            warnings.warn("EPW kind ('{0:}') not yet implemented - or it's unknown!".format(kind), Warning)

        self.pedestrian_wind_speed = pd.Series(index=self.index, name="pedestrian_wind_speed", data=wind_speed_at_height(source_wind_speed=self.wind_speed, source_wind_height=10, target_wind_height=1.5, terrain_roughness="Airport runway areas", log_method=True))

        try:
            weatherfile_ground_temperatures(self)
        except AttributeError as error:
            warnings.warn("{0:}\nGround temperatures not found in EPW file header".format(error), Warning)

    return self


def to_df(self):
    df = pd.concat([getattr(self, name) for name in dir(self) if type(getattr(self, name)).__name__ == "Series"], axis=1)
    return df


def to_wea(self, wea_file_path=None):
    """
    Create WEA file from Weather object.

    Parameters
    ----------
    self : Weather
        Weather object containing direct_normal_radiation and diffuse_horizontal_radiation
    wea_file_path : str
        The file path into which the WEA formatted radiation values will be written. If no value is passed, the
        output will be written to the same directory as the Weather object file_path.

    Returns
    -------
    str
        Path to the serialised WEA file.

    """

    if (self.direct_normal_radiation is None) | (self.diffuse_horizontal_radiation is None):
        raise Exception('No radiation data is available, try loading some first!')

    if wea_file_path is None:
        wea_file_path = pathlib.Path(self.file_path).with_suffix(".wea")

    header = "place {0:}_{1:}\nlatitude {2:0.3f}\nlongitude {3:0.3f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.1f}\nweather_data_file_units 1".format(
        slugify(self.city), slugify(self.country), self.latitude, -self.longitude, -self.time_zone * 15,
        self.elevation)

    values = []
    for n, dt in enumerate(self.index):
        values.append(
            "{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n],
                                                self.diffuse_horizontal_radiation[n]))
    with open(wea_file_path, "w") as f:
        f.write(header + "\n" + "\n".join(values) + "\n")

    self.wea_file = wea_file_path

    print("WEA file created: {0:}".format(wea_file_path.relative_to(wea_file_path.parent)))

    return self


def to_csv(self, csv_file_path=None):
    """
    Write DataFrame containing hourly annual weather variables, and derived variables to a CSV file.

    Parameters
    ----------
    csv_file_path : str
        The file path into which the annual hourly values will be written. If no value is passed, the output
        will be written to the same directory as the input weather-file.

    Returns
    -------
    str
        Path to the serialised CSV file.

    """

    df = self.to_df()

    if csv_file_path is None:
        csv_file_path = pathlib.Path(self.file_path).with_suffix(".csv")

    df.to_csv(csv_file_path)

    self.csv_file = csv_file_path

    print("CSV file created: {0:}".format(self.csv_file))

    return self


def to_pickle(self, pickle_path=None):

    if pickle_path is None:
        pickle_path = pathlib.Path(self.file_path).with_suffix(".pkl")
    self.pickle_path = pickle_path

    # Write pickle
    pickle.dump(self, open(self.pickle_path, "wb"))

    print("Weather object pickled: {0:}".format(self.pickle_path))

    return self


def load_pickle(pickle_path=None):
    return pickle.load(open(pickle_path, "rb"))
