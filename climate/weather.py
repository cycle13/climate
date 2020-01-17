import pathlib
import pickle
from io import StringIO
import pandas as pd

from .constants import DATETIME_INDEX, REINHART_PATCH_CENTROIDS, TREGENZA_PATCH_CENTROIDS, REINHART_PATCH_VECTORS, \
    REINHART_PATCH_COUNT, TREGENZA_PATCH_VECTORS, TREGENZA_PATCH_COUNT
from .helpers import slugify, chunk
from .psychrometrics import *
from .wind import *
from .sun import sun_position, create_sky_matrices
from .ground import interpolate_epw_ground_temperatures
from .material import OpaqueMaterial


def read_pickle(pickle_path=None):
    return pickle.load(open(pickle_path, "rb"))


class Weather(object):

    def __init__(self, epw_path=None):

        self.file_path = pathlib.Path(epw_path).absolute() if epw_path is not None else None
        self.wea_path = None
        self.csv_path = None
        self.pickle_path = None

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
        self.index = None

    def read_epw(self, kind="standard"):
        """Import annual hourly data from an EPW file

        Parameters
        ----------
        kind : str
            The kind of EPW being read. ["standard", "ccwwg", "cibse"]

        Returns
        -------
        Weather
        """

        assert self.file_path.is_file(), \
            'Cannot find EPW file at {}'.format(self.file_path)

        # Define the columns present in each different EPW kind
        epw_content = {
            "standard": ['year', 'month', 'day', 'hour', 'minute', 'data_source_and_uncertainty_flags',
                         'dry_bulb_temperature', 'dew_point_temperature', 'relative_humidity',
                         'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
                         'extraterrestrial_direct_normal_radiation', 'horizontal_infrared_radiation_intensity',
                         'global_horizontal_radiation', 'direct_normal_radiation', 'diffuse_horizontal_radiation',
                         'global_horizontal_illuminance', 'direct_normal_illuminance', 'diffuse_horizontal_illuminance',
                         'zenith_luminance', 'wind_direction', 'wind_speed', 'total_sky_cover', 'opaque_sky_cover',
                         'visibility', 'ceiling_height', 'present_weather_observation', 'present_weather_codes',
                         'precipitable_water', 'aerosol_optical_depth', 'snow_depth', 'days_since_last_snowfall',
                         'albedo', 'liquid_precipitation_depth', 'liquid_precipitation_quantity', ],
            "ccwwg": ['year', 'month', 'day', 'hour', 'minute', 'data_source_and_uncertainty_flags',
                      'dry_bulb_temperature', 'dew_point_temperature', 'relative_humidity',
                      'atmospheric_station_pressure', 'extraterrestrial_horizontal_radiation',
                      'extraterrestrial_direct_normal_radiation', 'horizontal_infrared_radiation_intensity',
                      'global_horizontal_radiation', 'direct_normal_radiation', 'diffuse_horizontal_radiation',
                      'global_horizontal_illuminance', 'direct_normal_illuminance', 'diffuse_horizontal_illuminance',
                      'zenith_luminance', 'wind_direction', 'wind_speed', 'total_sky_cover', 'opaque_sky_cover',
                      'visibility', 'ceiling_height', 'present_weather_observation', 'present_weather_codes',
                      'precipitable_water', 'aerosol_optical_depth', 'snow_depth', 'days_since_last_snowfall', ],
            "cibse": []
        }
        assert kind in epw_content.keys(), \
            "EPW kind '{}' has not been implemented".format(kind)

        with open(self.file_path, "r") as f:
            data = f.readlines()

            # Read location data
            self.city, self.region, self.country, self.dataset_type, self.station_id, self.latitude, self.longitude, self.time_zone, self.elevation = data[0].strip().split(",")[1:]
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

            # Read the data table and assign
            df = pd.read_csv(StringIO("\n".join(data[8:])), header=None).set_index(DATETIME_INDEX.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone))
            df.columns = epw_content[kind]
            df.drop(columns=["year", "month", "day", "hour", "minute"], inplace=True)

            # Set object attributes from DataFrame
            # TODO - This method is slow, but compact. A faster solution would be to assign to each attribute name manually, but then this makes the code look messy + needs handling of different EPW kinds.
            [setattr(self, col, df[col]) for col in df]
            self.index = df.index

            print("EPW loaded: {}".format(self.file_path))

        # The following lines automatically run to generate downstream variables
        self.calculate_psychrometrics()
        self.interpolate_ground_temperature()
        self.calculate_pedestrian_wind_speed()
        self.generate_sky_matrix()

        #
        # try:
        #     # Run psychrometrics on loaded EPW data
        #     self.calculate_psychrometrics()
        # except Exception as error:
        #     print("Looks like there isn't enough data to derive psychrometrics\n{}".format(error))

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

        self.wea_path = pathlib.Path(self.file_path).with_suffix(".wea") if wea_path is None else wea_path

        header = "place {0:}_{1:}\nlatitude {2:0.3f}\nlongitude {3:0.3f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.1f}\nweather_data_file_units 1".format(
            slugify(self.city), slugify(self.country), self.latitude, -self.longitude, -self.time_zone * 15,
            self.elevation)

        values = ["{0:} {1:} {2:}.5 {3:} {4:}".format(dt.month, dt.day, dt.hour, self.direct_normal_radiation.values[n], self.diffuse_horizontal_radiation[n]) for n, dt in enumerate(self.index)]

        with open(self.wea_path, "w") as f:
            f.write(header + "\n" + "\n".join(values) + "\n")

        print("WEA file created: {0:}".format(self.wea_path))

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

        self.csv_path = pathlib.Path(self.file_path).with_suffix(".csv") if csv_path is None else csv_path
        self.to_df().to_csv(self.csv_path)
        print("CSV file created: {0:}".format(self.csv_path))

    def to_pickle(self, pickle_path=None):
        self.pickle_path = pathlib.Path(self.file_path).with_suffix(".pkl") if pickle_path is None else pickle_path
        pickle.dump(self, open(self.pickle_path, "wb"))
        print("Object pickled: {0:}".format(self.pickle_path))

    def interpolate_ground_temperature(self):
        df = interpolate_epw_ground_temperatures(self.ground_temperatures, self.index)
        [setattr(self, col, df[col]) for col in df]
        return self

    def calculate_pedestrian_wind_speed(self, source_wind_height=10, target_wind_height=1.5, terrain_roughness="Airport runway areas", log_method=True):
        self.pedestrian_wind_speed = pd.Series(name="pedestrian_wind_speed", index=self.index, data=wind_speed_at_height(self.wind_speed, source_wind_height, target_wind_height, terrain_roughness, log_method))
        return self

    def calculate_psychrometrics(self):
        self.humidity_ratio = pd.Series(index=self.index, data=humidity_ratio_from_relative_humidity(self.dry_bulb_temperature, self.relative_humidity / 100, self.atmospheric_station_pressure))
        self.wet_bulb_temperature = pd.Series(index=self.index, data=wet_bulb_temperature_from_humidity_ratio(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure))
        self.partial_vapour_pressure_moist_air = pd.Series(index=self.index, data=vapor_pressure_from_humidity_ratio(self.humidity_ratio, self.atmospheric_station_pressure))
        self.enthalpy = pd.Series(index=self.index, data=enthalpy_from_moist_air(self.dry_bulb_temperature, self.humidity_ratio))
        self.specific_volume_moist_air = pd.Series(index=self.index, data=specific_volume_from_moist_air(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure))
        self.degree_of_saturation = pd.Series(index=self.index, data=degree_of_saturation(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure))

        return self

    def calculate_sun_position(self):

        df = sun_position(self.index, self.latitude, self.longitude)
        [setattr(self, col, df[col]) for col in df]
        return self

    def generate_sky_matrix(self, reinhart=True):

        # Create WEA file
        self.to_wea()

        # Calculate sky matrices
        self.direct_sky_matrix, self.diffuse_sky_matrix, self.total_sky_matrix = create_sky_matrices(self.wea_path, reinhart=reinhart)

        # Set the other descriptors for the sky patches
        self.sky_matrix_patch_centroids = REINHART_PATCH_CENTROIDS if reinhart else TREGENZA_PATCH_CENTROIDS
        self.sky_matrix_patch_vectors = REINHART_PATCH_VECTORS if reinhart else TREGENZA_PATCH_VECTORS
        self.sky_matrix_patch_count = REINHART_PATCH_COUNT if reinhart else TREGENZA_PATCH_COUNT

        return self

    def calculate_mrt(self, method="openfield", ground_material=None):

        # if ground_material is None:
        #     ground_material = OpaqueMaterial(specific_heat_capacity=, reflectivity=, emissivity=, thickness=, roughness="Medium rough")
        return self

