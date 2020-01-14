from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, \
    GetMoistAirEnthalpy, GetHumRatioFromEnthalpyAndTDryBulb, GetMoistAirVolume, GetDegreeOfSaturation, GetHumRatioFromVapPres, GetTDewPointFromVapPres, \
    GetTDryBulbFromEnthalpyAndHumRatio, GetSatVapPres, GetTDryBulbFromMoistAirVolumeAndHumRatio, GetDryAirEnthalpy, GetSatAirEnthalpy

import pandas as pd
import numpy as np

SetUnitSystem(SI)

def humidity_ratio_from_partial_vapor_pressure(water_partial_vapor_pressure, pressure):
    f_GetHumRatioFromVapPres = np.vectorize(GetHumRatioFromVapPres)
    return f_GetHumRatioFromVapPres(water_partial_vapor_pressure, pressure)


def humidity_ratio_from_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
    f_GetHumRatioFromRelHum = np.vectorize(GetHumRatioFromRelHum)
    return f_GetHumRatioFromRelHum(dry_bulb_temperature, relative_humidity / 100, pressure)


def humidity_ratio_from_enthalpy_moist_air(moist_air_enthalpy, dry_bulb_temperature):
    f_GetHumRatioFromEnthalpyAndTDryBulb = np.vectorize(GetHumRatioFromEnthalpyAndTDryBulb)
    return f_GetHumRatioFromEnthalpyAndTDryBulb(moist_air_enthalpy, dry_bulb_temperature)


def enthalpy_from_dry_air(dry_bulb_temperature):
    f_GetDryAirEnthalpy = np.vectorize(GetDryAirEnthalpy)
    return f_GetDryAirEnthalpy(dry_bulb_temperature)


def enthalpy_from_moist_air(dry_bulb_temperature, humidity_ratio):
    f_GetMoistAirEnthalpy = np.vectorize(GetMoistAirEnthalpy)
    return f_GetMoistAirEnthalpy(dry_bulb_temperature, humidity_ratio)


def enthalpy_from_saturated_air(dry_bulb_temperature, pressure):
    f_GetSatAirEnthalpy = np.vectorize(GetSatAirEnthalpy)
    return f_GetSatAirEnthalpy(dry_bulb_temperature, pressure)


def specific_volume_from_moist_air(dry_bulb_temperature, humidity_ratio, pressure):
    f_GetMoistAirVolume = np.vectorize(GetMoistAirVolume)
    return f_GetMoistAirVolume(dry_bulb_temperature, humidity_ratio, pressure)


def dew_point_temperature_from_vapor_pressue(dry_bulb_temperature, vapor_pressure):
    f_GetTDewPointFromVapPres = np.vectorize(GetTDewPointFromVapPres)
    return f_GetTDewPointFromVapPres(dry_bulb_temperature, vapor_pressure)


def dry_bulb_temperature_from_enthalpy(moist_air_enthalpy, humidity_ratio):
    f_GetTDryBulbFromEnthalpyAndHumRatio = np.vectorize(GetTDryBulbFromEnthalpyAndHumRatio)
    return f_GetTDryBulbFromEnthalpyAndHumRatio(moist_air_enthalpy, humidity_ratio)


def wet_bulb_temperature_from_relative_humidity(dry_bulb_temperature, relative_humidity, pressure):
    f_GetTWetBulbFromRelHum = np.vectorize(GetTWetBulbFromRelHum)
    return f_GetTWetBulbFromRelHum(dry_bulb_temperature, relative_humidity, pressure)


def wet_bulb_temperature_from_humidity_ratio(dry_bulb_temperature, humidity_ratio, pressure):
    f_GetTWetBulbFromHumRatio = np.vectorize(GetTWetBulbFromHumRatio)
    return f_GetTWetBulbFromHumRatio(dry_bulb_temperature, humidity_ratio, pressure)


def saturation_vapor_pressure(dry_bulb_temperature):
    f_GetSatVapPres = np.vectorize(GetSatVapPres)
    return f_GetSatVapPres(dry_bulb_temperature)


def degree_of_saturation(dry_bulb_temperature, humidity_ratio, pressure):
    f_GetDegreeOfSaturation = np.vectorize(GetDegreeOfSaturation)
    return f_GetDegreeOfSaturation(dry_bulb_temperature, humidity_ratio, pressure)


def vapor_pressure_from_humidity_ratio(humidity_ratio, pressure):
    f_GetVapPresFromHumRatio = np.vectorize(GetVapPresFromHumRatio)
    return f_GetVapPresFromHumRatio(humidity_ratio, pressure)


def dry_bulb_temperature_from_moist_air_volume_and_humidity_ratio(moist_air_volume, humidity_ratio_pressure):
    f_GetTDryBulbFromMoistAirVolumeAndHumRatio = np.vectorize(GetTDryBulbFromMoistAirVolumeAndHumRatio)
    return f_GetTDryBulbFromMoistAirVolumeAndHumRatio()


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
    enthalpy = enthalpy_from_moist_air(dry_bulb_temperature, humidity_ratio)
    specific_volume_moist_air = specific_volume_from_moist_air(dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure)
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

