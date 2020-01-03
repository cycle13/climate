from psychrochart.psychrolib_extra import GetTDryBulbFromMoistAirVolume
from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, \
    GetMoistAirEnthalpy, GetHumRatioFromEnthalpyAndTDryBulb, GetMoistAirVolume, GetDegreeOfSaturation, GetHumRatioFromVapPres, GetTDewPointFromVapPres, \
    GetTDryBulbFromEnthalpyAndHumRatio, GetSatVapPres, GetTDryBulbFromMoistAirVolumeAndHumRatio, GetDryAirEnthalpy, GetSatAirEnthalpy
import pandas as pd
import numpy as np

SetUnitSystem(SI)

# Matricize computation methods
humidity_ratio_from_vapor_pressure = np.vectorize(GetHumRatioFromVapPres)
humidity_ratio_from_relative_humidity = np.vectorize(GetHumRatioFromRelHum)
humidity_ratio_from_enthalpy = np.vectorize(GetHumRatioFromEnthalpyAndTDryBulb)
dry_air_enthalpy = np.vectorize(GetDryAirEnthalpy)
moist_air_enthalpy = np.vectorize(GetMoistAirEnthalpy)
sat_air_enthalpy = np.vectorize(GetSatAirEnthalpy)
moist_air_volume = np.vectorize(GetMoistAirVolume)
dew_point_from_vapor_pressue = np.vectorize(GetTDewPointFromVapPres)
dry_bulb_temperature_from_enthalpy = np.vectorize(GetTDryBulbFromEnthalpyAndHumRatio)
dry_bulb_temperature_from_specific_volume = np.vectorize(GetTDryBulbFromMoistAirVolume)
wet_bulb_temperature_from_humidity_ratio = np.vectorize(GetTWetBulbFromHumRatio)
saturation_vapor_pressure = np.vectorize(GetSatVapPres)
degree_of_saturation = np.vectorize(GetDegreeOfSaturation)
vapor_pressure_from_humidity_ratio = np.vectorize(GetVapPresFromHumRatio)
dry_bulb_temperature_from_moist_air_volume_and_humidity_ratio = np.vectorize(GetTDryBulbFromMoistAirVolumeAndHumRatio)

def psychrometric_calculations(dry_bulb_temperature, relative_humidity, atmospheric_station_pressure):
    """
    Calculate a range of derivable air/water characteristics from Dry-bulb temperature, Relative humidity and Atmospheric pressure.

    Parameters
    ----------
    dry_bulb_temperature : float
        Dry-bulb temperature of air in degrees C
    relative_humidity : float
        Relative humidity of air in %
    atmospheric_station_pressure : float
        Atmospheric pressure in Pa

    Returns
    -------
    DataFrame
        Pandas DataFrame object containing calculated psychrometric values

    """

    SetUnitSystem(SI)

    # Check for inputs of single floats
    if type(dry_bulb_temperature) == type(relative_humidity) == type(atmospheric_station_pressure) == int:
        dry_bulb_temperature = np.array([dry_bulb_temperature])
        relative_humidity = np.array([relative_humidity])
        atmospheric_station_pressure = np.array([atmospheric_station_pressure])

    # Convert inputs into numpy arrays for efficient downstream handling
    dry_bulb_temperature = np.array(dry_bulb_temperature)
    relative_humidity = np.array(relative_humidity)
    atmospheric_station_pressure = np.array(atmospheric_station_pressure)

    # Check that inputs are all the same length
    if not len(dry_bulb_temperature) == len(relative_humidity) == len(atmospheric_station_pressure):
        raise ValueError("Inputs should all be the same length")

    # Calculate psychrometrics
    humidity_ratio = humidity_ratio_from_relative_humidity(dry_bulb_temperature, relative_humidity / 100, atmospheric_station_pressure)
    wet_bulb_temperature = wet_bulb_temperature_from_humidity_ratio(dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure)
    partial_vapour_pressure_moist_air = vapor_pressure_from_humidity_ratio(humidity_ratio, atmospheric_station_pressure)
    enthalpy = moist_air_enthalpy(dry_bulb_temperature, humidity_ratio)
    specific_volume_moist_air = moist_air_volume(dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure)
    deg_saturation = degree_of_saturation(dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure)

    return humidity_ratio, wet_bulb_temperature, partial_vapour_pressure_moist_air, enthalpy, specific_volume_moist_air, deg_saturation


def annual_psychrometrics(self):
    humidity_ratio, wet_bulb_temperature, partial_vapour_pressure_moist_air, enthalpy, specific_volume_moist_air, degree_of_saturation = psychrometric_calculations(self.dry_bulb_temperature, self.relative_humidity, self.atmospheric_station_pressure)
    self.humidity_ratio = pd.Series(name="humidity_ratio", index=self.index, data=humidity_ratio)
    self.wet_bulb_temperature = pd.Series(name="wet_bulb_temperature", index=self.index, data=wet_bulb_temperature)
    self.partial_vapour_pressure_moist_air = pd.Series(name="partial_vapour_pressure_moist_air", index=self.index, data=partial_vapour_pressure_moist_air)
    self.enthalpy = pd.Series(name="enthalpy", index=self.index, data=enthalpy)
    self.specific_volume_moist_air = pd.Series(name="specific_volume_moist_air", index=self.index, data=specific_volume_moist_air)
    self.degree_of_saturation = pd.Series(name="degree_of_saturation", index=self.index, data=degree_of_saturation)
    print("Psychrometric calculations successful")
    return self.humidity_ratio, self.wet_bulb_temperature, self.partial_vapour_pressure_moist_air, self.enthalpy, self.specific_volume_moist_air, self.degree_of_saturation

