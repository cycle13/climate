import pathlib
import pandas as pd
import numpy as np
from ladybug.epw import EPW

from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, GetMoistAirEnthalpy, GetMoistAirVolume, GetDegreeOfSaturation


class Weather(object):

    def __init__(self, epw_path, sky_matrix=None):

        # Instantiate the objects core properties
        self.epw_path = pathlib.Path(epw_path).absolute() if epw_path is not None else None

        # Load weatherfile
        epw = EPW(self.epw_path)

        # Assign header information to variables
        self.city = epw.location.city
        self.country = epw.location.country
        self.elevation = epw.location.elevation
        self.latitude = epw.location.latitude
        self.longitude = epw.location.longitude
        self.meridian = epw.location.meridian
        self.source = epw.location.source
        self.state = epw.location.state
        self.station_id = epw.location.station_id
        self.time_zone = epw.location.time_zone

        # Define datetime index based on UTC-offset
        idx_month = pd.date_range(start="2018-01-01 00:00:00", end="2019-01-01 00:00:00", freq="1M", closed="left").tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone)
        idx_year = pd.date_range(start="2018-01-01 00:00:00", end="2019-01-01 00:00:00", freq="60T", closed="left").tz_localize("UTC").tz_convert(int(self.time_zone * 60 * 60)) - pd.Timedelta(hours=self.time_zone)

        # Ground temperature data
        for depth in [0.5, 2, 4]:
            gt = pd.Series(index=idx_month, data=epw.monthly_ground_temperature[depth].values, name="ground_temperature_{0:}mm".format(depth * 1000)).resample("MS").mean().reindex(idx_year)
            gt.iloc[-1] = gt.iloc[0]
            gt.interpolate(inplace=True)
            setattr(self, "ground_temperature_{0:}mm".format(depth * 1000), gt)

        # Add climate variables
        self.dry_bulb_temperature = pd.Series(index=idx_year, data=epw.dry_bulb_temperature.values, name="dry_bulb_temperature")
        self.dew_point_temperature = pd.Series(index=idx_year, data=epw.dew_point_temperature.values, name="dew_point_temperature")
        self.relative_humidity = pd.Series(index=idx_year, data=epw.relative_humidity.values, name="relative_humidity")
        self.atmospheric_station_pressure = pd.Series(index=idx_year, data=epw.atmospheric_station_pressure.values, name="atmospheric_station_pressure")
        self.extraterrestrial_horizontal_radiation = pd.Series(index=idx_year, data=epw.extraterrestrial_horizontal_radiation.values, name="extraterrestrial_horizontal_radiation")
        self.extraterrestrial_direct_normal_radiation = pd.Series(index=idx_year, data=epw.extraterrestrial_direct_normal_radiation.values, name="extraterrestrial_direct_normal_radiation")
        self.horizontal_infrared_radiation_intensity = pd.Series(index=idx_year, data=epw.horizontal_infrared_radiation_intensity.values, name="horizontal_infrared_radiation_intensity")
        self.global_horizontal_radiation = pd.Series(index=idx_year, data=epw.global_horizontal_radiation.values, name="global_horizontal_radiation")
        self.direct_normal_radiation = pd.Series(index=idx_year, data=epw.direct_normal_radiation.values, name="direct_normal_radiation")
        self.diffuse_horizontal_radiation = pd.Series(index=idx_year, data=epw.diffuse_horizontal_radiation.values, name="diffuse_horizontal_radiation")
        self.global_horizontal_illuminance = pd.Series(index=idx_year, data=epw.global_horizontal_illuminance.values, name="global_horizontal_illuminance")
        self.direct_normal_illuminance = pd.Series(index=idx_year, data=epw.direct_normal_illuminance.values, name="direct_normal_illuminance")
        self.diffuse_horizontal_illuminance = pd.Series(index=idx_year, data=epw.diffuse_horizontal_illuminance.values, name="diffuse_horizontal_illuminance")
        self.zenith_luminance = pd.Series(index=idx_year, data=epw.zenith_luminance.values, name="zenith_luminance")
        self.wind_direction = pd.Series(index=idx_year, data=epw.wind_direction.values, name="wind_direction")
        self.wind_speed = pd.Series(index=idx_year, data=epw.wind_speed.values, name="wind_speed")
        self.total_sky_cover = pd.Series(index=idx_year, data=epw.total_sky_cover.values, name="total_sky_cover")
        self.opaque_sky_cover = pd.Series(index=idx_year, data=epw.opaque_sky_cover.values, name="opaque_sky_cover")
        self.visibility = pd.Series(index=idx_year, data=epw.visibility.values, name="visibility")
        self.ceiling_height = pd.Series(index=idx_year, data=epw.ceiling_height.values, name="ceiling_height")
        self.present_weather_observation = pd.Series(index=idx_year, data=epw.present_weather_observation.values, name="present_weather_observation")
        self.present_weather_codes = pd.Series(index=idx_year, data=epw.present_weather_codes.values, name="present_weather_codes")
        self.precipitable_water = pd.Series(index=idx_year, data=epw.precipitable_water.values, name="precipitable_water")
        self.aerosol_optical_depth = pd.Series(index=idx_year, data=epw.aerosol_optical_depth.values, name="aerosol_optical_depth")
        self.snow_depth = pd.Series(index=idx_year, data=epw.snow_depth.values, name="snow_depth")
        self.days_since_last_snowfall = pd.Series(index=idx_year, data=epw.days_since_last_snowfall.values, name="days_since_last_snowfall")
        self.albedo = pd.Series(index=idx_year, data=epw.albedo.values, name="albedo")
        self.liquid_precipitation_depth = pd.Series(index=idx_year, data=epw.liquid_precipitation_depth.values, name="liquid_precipitation_depth")
        self.liquid_precipitation_quantity = pd.Series(index=idx_year, data=epw.liquid_precipitation_quantity.values, name="liquid_precipitation_quantity")
        self.sky_temperature = pd.Series(index=idx_year, data=epw.sky_temperature.values, name="sky_temperature")

        # Add derived climate variables (using Psychrolib for speed)
        SetUnitSystem(SI)
        self.humidity_ratio = pd.Series(index=idx_year, data=np.vectorize(GetHumRatioFromRelHum)(self.dry_bulb_temperature, self.relative_humidity / 100, self.atmospheric_station_pressure), name="humidity_ratio")
        self.wet_bulb_temperature = pd.Series(index=idx_year, data=np.vectorize(GetTWetBulbFromHumRatio)(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure), name="wet_bulb_temperature")
        self.partial_vapour_pressure_moist_air = pd.Series(index=idx_year, data=np.vectorize(GetVapPresFromHumRatio)(self.humidity_ratio, self.atmospheric_station_pressure), name="partial_vapour_pressure_moist_air")
        self.enthalpy = pd.Series(index=idx_year, data=np.vectorize(GetMoistAirEnthalpy)(self.dry_bulb_temperature, self.humidity_ratio), name="enthalpy")
        self.specific_volume_moist_air = pd.Series(index=idx_year, data=np.vectorize(GetMoistAirVolume)(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure), name="specific_volume_moist_air")
        self.degree_of_saturation = pd.Series(index=idx_year, data=np.vectorize(GetDegreeOfSaturation)(self.dry_bulb_temperature, self.humidity_ratio, self.atmospheric_station_pressure), name="degree_of_saturation")
