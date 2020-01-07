from climate.common.helpers import *
from climate.compute.psychrometrics import annual_psychrometrics
from climate.compute.sun import annual_sun_position, generate_sky_matrix
from climate.compute.mean_radiant_temperature import mrt_solar_adjusted, mrt_openfield
from climate.compute.comfort import utci_openfield, \
    utci_solar_adjusted, set_solar_adjusted, set_openfield
from climate.compute.ground_temperature import weatherfile_ground_temperatures, annual_ground_temperature_at_depth
from climate.compute.wind import pedestrian_wind_speed
from climate.compute.uwg import run_uwg

from climate.plot.diurnal import diurnal
from climate.plot.heatmap import heatmap
from climate.plot.radiation_rose import radiation_rose
from climate.plot.wind_rose import windrose
from climate.plot.psychrometric import psychrometric
from climate.plot.utci import utci_frequency, utci_heatmap

import pandas as pd
import pathlib
from io import StringIO


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
        self.index = None
        # self.df = None
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

        # Loaded and calculated ground temperatures
        self.ground_temperature_500_weatherfile = None
        self.ground_temperature_2000_weatherfile = None
        self.ground_temperature_4000_weatherfile = None
        self.ground_temperature_calculated = None

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
        self.mean_radiant_temperature_solar_adjusted = None
        self.mean_radiant_temperature_openfield = None

        self.universal_thermal_climate_index_solar_adjusted = None
        self.universal_thermal_climate_index_openfield = None

        self.standard_effective_temperature_solar_adjusted = None
        self.standard_effective_temperature_openfield = None

        # Urban weather generator values
        self.uwg_parameters = None
        self.uwg_weatherfile = None

    def read(self, sun_position=True, pedestrian_wind=True, psychrometrics=True, sky_matrix=False, ground_temp=False, mrt=False, utci=False, _set=False, uwg=False):
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
            self.index = DATETIME_INDEX.tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone)
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
            # self.df = df

        # Generate the sky matrix
        if sky_matrix:
            generate_sky_matrix(self)

        if ground_temp:
            # Calculate soil temperatures
            annual_ground_temperature_at_depth(self, depth=0.5)

            # Interpolate ground temperatures from loaded weatherfile
            weatherfile_ground_temperatures(self)

        if pedestrian_wind:
            # Calculate pedestrian height wind speed (at 1.5m above ground)
            pedestrian_wind_speed(self)

        if sun_position:
            # Run the solar position calculations
            annual_sun_position(self)

        if psychrometrics:
            # Run the psychrometric calculations
            annual_psychrometrics(self)

        if mrt:
            # Run the Solar adjusted MRT method
            mrt_solar_adjusted(self)

            # Run the openfield MRT method
            mrt_openfield(self)

        if utci:

            # Run the UTCI using the solar adjusted MRT values
            utci_solar_adjusted(self)

            # Run the UTCI using the openfield MRT values
            utci_openfield(self)

        if _set:
            # Run the SET using the solar adjusted MRT values
            set_solar_adjusted(self)

            # Run the SET using the openfield MRT values
            set_openfield(self)

        if uwg:
            # Run the urban weather generator!
            run_uwg(self)

        return self

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
        header = "place {0:}_{1:}\nlatitude {2:0.3f}\nlongitude {3:0.3f}\ntime_zone {4:0.2f}\nsite_elevation {5:0.1f}\nweather_data_file_units 1".format(
            slugify(self.city), slugify(self.country), self.latitude, -self.longitude, -self.time_zone * 15,
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

    def to_df(self):
        df = pd.concat([getattr(self, name) for name in dir(self) if type(getattr(self, name)).__name__ == "Series"], axis=1)
        return df

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
        df = self.to_df()
        if file_path is None:
            file_path = pathlib.Path(self.file_path).with_suffix(".csv")
        df.to_csv(file_path)
        self.csv_file = file_path
        print("CSV file created: {0:}".format(self.csv_file))
        return df

    def plot(self):

        diurnal(self, dew_point=True, tone_color="#555555", save=True)
        diurnal(self, dew_point=False, tone_color="#555555", save=True)

        heatmap(self, "dry_bulb_temperature", cmap="Reds", tone_color="#555555", save=True)
        heatmap(self, "relative_humidity", cmap="Blues", tone_color="#555555", save=True)
        heatmap(self, "wind_direction", cmap="Greys", tone_color="#555555", save=True)
        heatmap(self, "direct_normal_radiation", cmap="Oranges", tone_color="#555555", save=True)
        heatmap(self, "diffuse_horizontal_radiation", cmap="Oranges", tone_color="#555555", save=True)
        heatmap(self, "global_horizontal_radiation", cmap="Oranges", tone_color="#555555", save=True)

        psychrometric(self, bins=50, cmap="inferno", tone_color="#555555", save=True)

        for season_period in ["Annual"]:#, "Spring", "Summer", "Autumn", "Winter"]:
            for day_period in ["Daily"]:#, "Morning", "Midday", "Afternoon", "Evening", "Night"]:

                radiation_rose(self, season_period=season_period, day_period=day_period, n_sector=36, cmap=None, tone_color="#555555", same_scale=False, save=True)
                windrose(self, season_period=season_period, day_period=day_period, n_sector=16, cmap=None, tone_color="#555555", save=True)

        try:
            utci_frequency(self, "universal_thermal_climate_index_openfield", tone_color="#555555", save=True)
            utci_frequency(self, "universal_thermal_climate_index_solar_adjusted", tone_color="#555555", save=True)

            utci_heatmap(self, "universal_thermal_climate_index_openfield", tone_color="#555555", save=True)
            utci_heatmap(self, "universal_thermal_climate_index_solar_adjusted", tone_color="#555555", save=True)
        except Exception as e:
            print("UTCI hasn't been run - or something like that. I don't know. This bit hasn't been sorted out yet\n{0:}".format(e))



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
