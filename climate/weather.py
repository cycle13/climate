import pathlib
import pandas as pd
import numpy as np
from io import StringIO
from pvlib.solarposition import get_solarposition

from .helpers import slugify, chunk

from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, GetMoistAirEnthalpy, GetMoistAirVolume, GetDegreeOfSaturation
SetUnitSystem(SI)


class Weather(object):

    def __init__(self, epw_path=None):

        self.epw_path = pathlib.Path(epw_path).absolute() if epw_path is not None else None
        self.index = pd.date_range(start="2018-01-01 00:00:00", end="2019-01-01 00:00:00", freq="60T", closed="left")
        self.read_epw()

    def read_epw(self):
        """ Import annual hourly data from an EPW file
        """

        assert self.epw_path.is_file(), 'Cannot find EPW file at {}'.format(self.epw_path)

        epw_content = ['year', 'month', 'day', 'hour', 'minute', 'data_source_and_uncertainty_flags',
                         'dry_bulb_temperature', 'dew_point_temperature', 'relative_humidity',
                         'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
                         'extraterrestrial_direct_normal_radiation', 'horizontal_infrared_radiation_intensity',
                         'global_horizontal_radiation', 'direct_normal_radiation', 'diffuse_horizontal_radiation',
                         'global_horizontal_illuminance', 'direct_normal_illuminance', 'diffuse_horizontal_illuminance',
                         'zenith_luminance', 'wind_direction', 'wind_speed', 'total_sky_cover', 'opaque_sky_cover',
                         'visibility', 'ceiling_height', 'present_weather_observation', 'present_weather_codes',
                         'precipitable_water', 'aerosol_optical_depth', 'snow_depth', 'days_since_last_snowfall',
                         'albedo', 'liquid_precipitation_depth', 'liquid_precipitation_quantity', ]

        with open(self.epw_path, "r") as f:
            data = f.readlines()

            # Read location data
            city, region, country, dataset_type, station_id, latitude, longitude, time_zone, elevation = data[0].strip().split(",")[1:]
            setattr(self, "city", city)
            setattr(self, "region", region)
            setattr(self, "country", country)
            setattr(self, "dataset_type", dataset_type)
            setattr(self, "station_id", station_id)
            setattr(self, "latitude", float(latitude))
            setattr(self, "longitude", float(longitude))
            setattr(self, "time_zone", float(time_zone))
            setattr(self, "elevation", float(elevation))
            setattr(self, "design_conditions", ",".join(data[1].strip().split(",")[1:]))
            setattr(self, "typical_extreme_periods", ",".join(data[2].strip().split(",")[1:]))
            header_ground_temperatures = ",".join(data[3].strip().split(",")[1:])
            setattr(self, "ground_temperatures", header_ground_temperatures)
            setattr(self, "holidays_daylight_savings", ",".join(data[4].strip().split(",")[1:]))
            setattr(self, "comments_1", ",".join(data[5].strip().split(",")[1:]))
            setattr(self, "comments_2", ",".join(data[6].strip().split(",")[1:]))

            idx = self.index.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone)

            # Read the data table and assign
            df = pd.read_csv(StringIO("\n".join(data[8:])), header=None).set_index(idx)
            df.columns = epw_content
            df.drop(columns=["year", "month", "day", "hour", "minute", "data_source_and_uncertainty_flags"], inplace=True)

            # Add interpolated ground temperatures
            g_temps = {}
            for n, i in enumerate(list(chunk(header_ground_temperatures.split(",")[1:], n=16, method="size"))):
                g_temps[float(n)] = [float(j) for j in i[4:]]
            gtdf = pd.DataFrame.from_dict(g_temps)
            gtdf.index = pd.Series(index=idx).resample("MS").mean().index
            gtdf = pd.concat([pd.DataFrame(index=idx), gtdf], axis=1)
            gtdf.columns = ["ground_temperature_500mm", "ground_temperature_2000mm", "ground_temperature_4000mm"]
            gtdf.iloc[-1, :] = gtdf.iloc[0, :]  # Assign start temp to last datetime
            gtdf.interpolate(inplace=True)  # Fill in the gaps
            df = pd.concat([df, gtdf], axis=1)

            # Add sun position values
            spdf = get_solarposition(self.index, self.latitude, self.longitude)
            spdf.rename(
                columns={
                    'apparent_zenith': 'solar_apparent_zenith_angle',
                    'zenith': 'solar_zenith_angle',
                    'apparent_elevation': 'solar_apparent_elevation_angle',
                    'elevation': 'solar_elevation_angle',
                    'azimuth': 'solar_azimuth_angle',
                    'equation_of_time': 'solar_equation_of_time'
                }, inplace=True)
            spdf.index = idx
            df = pd.concat([df, spdf], axis=1)

            # Add psychrometrics
            df["humidity_ratio"] = np.vectorize(GetHumRatioFromRelHum)(df.dry_bulb_temperature, df.relative_humidity / 100, df.atmospheric_station_pressure)
            df["wet_bulb_temperature"] = np.vectorize(GetTWetBulbFromHumRatio)(df.dry_bulb_temperature, df.humidity_ratio, df.atmospheric_station_pressure)
            df["partial_vapour_pressure_moist_air"] = np.vectorize(GetVapPresFromHumRatio)(df.humidity_ratio, df.atmospheric_station_pressure)
            df["enthalpy"] = np.vectorize(GetMoistAirEnthalpy)(df.dry_bulb_temperature, df.humidity_ratio)
            df["specific_volume_moist_air"] = np.vectorize(GetMoistAirVolume)(df.dry_bulb_temperature, df.humidity_ratio, df.atmospheric_station_pressure)
            df["deg_saturation"] = np.vectorize(GetDegreeOfSaturation)(df.dry_bulb_temperature, df.humidity_ratio, df.atmospheric_station_pressure)

            # Get pedestrian wind speed
            df["pedestrian_wind_speed"] = wind_speed_at_height(df.wind_speed, 10, 1)

            # Set object attributes dynamically from DataFrame
            [setattr(self, col, df[col]) for col in df]

            print("EPW loaded: [{0:}]".format(self.epw_path))

        return self

    def to_df(self):
        return pd.concat([getattr(self, name) for name in dir(self) if type(getattr(self, name)).__name__ == "Series"], axis=1)

    def to_wea(self, wea_path=None):
        """ Create WEA file from Weather object.

        Parameters
        ----------
        self : Weather
            Weather object containing direct_normal_radiation and diffuse_horizontal_radiation
        wea_path : str
            The file path into which the WEA formatted radiation values will be written. If no value is passed, the
            output will be written to the same directory as the Weather object file_path.

        Returns
        -------
        str
            Path to the serialised WEA file.

        """

        if (self.direct_normal_radiation is None) | (self.diffuse_horizontal_radiation is None):
            raise Exception('No radiation data is available, try loading some first!')

        self.wea_path = pathlib.Path(self.epw_path).with_suffix(".wea") if wea_path is None else wea_path

        header = "place {0:}_{1:}\nlatitude {2:0.3f}\nlongitude {3:0.3f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.1f}\nweather_data_file_units 1".format(
            slugify(self.city), slugify(self.country), self.latitude, -self.longitude, -self.time_zone * 15,
            self.elevation)

        values = ["{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n], self.diffuse_horizontal_radiation[n]) for n, dt in enumerate(self.index)]

        with open(self.wea_path, "w") as f:
            f.write(header + "\n" + "\n".join(values) + "\n")

        print("WEA file created: {0:}".format(self.wea_path))
        return self.wea_path

    def to_csv(self, csv_path=None):
        """
        Write DataFrame containing hourly annual weather variables, and derived variables to a CSV file.

        Parameters
        ----------
        csv_path : str
            The file path into which the annual hourly values will be written. If no value is passed, the output
            will be written to the same directory as the input weather-file.

        Returns
        -------
        str
            Path to the serialised CSV file.

        """

        self.csv_path = pathlib.Path(self.epw_path).with_suffix(".csv") if csv_path is None else csv_path
        self.to_df().to_csv(self.csv_path)
        print("CSV file created: {0:}".format(self.csv_path))


def wind_speed_at_height(source_wind_speed, source_wind_height, target_wind_height, terrain_roughness="Airport runway areas", log_method=True):
    """

    :param source_wind_speed: The wind speed to be translated (m/s)
    :param source_wind_height: The height at which the source wind speed was measured
    :param target_wind_height: The height to which the source wind speed will be translated
    :param terrain_roughness: A terrain roughness value from the European Wind Atlas, page 58, Figure 3.1: Roughness length, surface characteristics and roughness class.
    :param log_method:
    :return:
    """
    roughness = {
        "City": 1,
        "Forest": 0.8,
        "Suburbs": 0.5,
        "Shelter belts": 0.3,
        "Many trees and/or bushes": 0.2,
        "Farmland with closed appearance": 0.1,
        "Farmland with open appearance": 0.05,
        "Farmland with very few buildings, trees etc. airport areas with buildings and trees": 0.03,
        "Airport runway areas": 0.01,
        "Mown grass": 0.0075,
        "Bare soil (smooth)": 0.005,
        "Snow surfaces (smooth)": 0.001,
        "Sand surfaces (smooth)": 0.0003,
        "Water areas (lakes, fjords, open sea)": 0.0001
    }

    if log_method:
        values = np.where(source_wind_speed <= 0, 0, source_wind_speed * (
                    np.log(target_wind_height / roughness[terrain_roughness]) / np.log(
                source_wind_height / roughness[terrain_roughness])))
    else:
        wind_shear_exponent = 1 / 7
        values = np.where(source_wind_speed <= 0, 0,
                          source_wind_speed * np.power(target_wind_height / source_wind_height, wind_shear_exponent))

    return values
