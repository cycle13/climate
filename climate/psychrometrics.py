from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, GetMoistAirEnthalpy, GetMoistAirVolume, GetDegreeOfSaturation
import pandas as pd
import numpy as np

SetUnitSystem(SI)

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

    # Wrap calculation methods to matricize computation
    def humidity_ratio_from_relative_humidity(input_variables):
        return GetHumRatioFromRelHum(input_variables[0], input_variables[1], input_variables[2])

    def wet_bulb_temperature_from_humidity_ratio(input_variables):
        return GetTWetBulbFromHumRatio(input_variables[0], input_variables[1], input_variables[2])

    def vapour_pressure_from_humidity_ratio(input_variables):
        return GetVapPresFromHumRatio(input_variables[0], input_variables[1])

    def enthalpy_from_humidity_ratio(input_variables):
        return GetMoistAirEnthalpy(input_variables[0], input_variables[1])

    def specific_volume_moist_air_from_humidity_ratio(input_variables):
        return GetMoistAirVolume(input_variables[0], input_variables[1], input_variables[2])

    def degree_of_saturation_from_humidity_ratio(input_variables):
        return GetDegreeOfSaturation(input_variables[0], input_variables[1], input_variables[2])

    # Calculate psychrometrics
    humidity_ratio = np.apply_along_axis(humidity_ratio_from_relative_humidity, 0, np.array([dry_bulb_temperature, relative_humidity / 100, atmospheric_station_pressure]))
    wet_bulb_temperature = np.apply_along_axis(wet_bulb_temperature_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))
    partial_vapour_pressure_moist_air = np.apply_along_axis(vapour_pressure_from_humidity_ratio, 0, np.array([humidity_ratio, atmospheric_station_pressure]))
    enthalpy = np.apply_along_axis(enthalpy_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio]))
    specific_volume_moist_air = np.apply_along_axis(specific_volume_moist_air_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))
    degree_of_saturation = np.apply_along_axis(degree_of_saturation_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))

    return humidity_ratio, wet_bulb_temperature, partial_vapour_pressure_moist_air, enthalpy, specific_volume_moist_air, degree_of_saturation


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

