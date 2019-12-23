from psychrolib import SetUnitSystem, SI, GetHumRatioFromRelHum, GetTWetBulbFromHumRatio, GetVapPresFromHumRatio, GetMoistAirEnthalpy, GetMoistAirVolume, GetDegreeOfSaturation
import pandas as pd
import numpy as np


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
    def humidity_ratio_from_relative_humidity(input):
        return GetHumRatioFromRelHum(input[0], input[1], input[2])

    def wet_bulb_temperature_from_humidity_ratio(input):
        return GetTWetBulbFromHumRatio(input[0], input[1], input[2])

    def vapour_pressure_from_humidity_ratio(input):
        return GetVapPresFromHumRatio(input[0], input[1])

    def enthalpy_from_humidity_ratio(input):
        return GetMoistAirEnthalpy(input[0], input[1])

    def specific_volume_moist_air_from_humidity_ratio(input):
        return GetMoistAirVolume(input[0], input[1], input[2])

    def degree_of_saturation_from_humidity_ratio(input):
        return GetDegreeOfSaturation(input[0], input[1], input[2])

    # Calculate psychrometrics
    humidity_ratio = np.apply_along_axis(humidity_ratio_from_relative_humidity, 0, np.array([dry_bulb_temperature, relative_humidity / 100, atmospheric_station_pressure]))
    wet_bulb_temperature = np.apply_along_axis(wet_bulb_temperature_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))
    partial_vapour_pressure_moist_air = np.apply_along_axis(vapour_pressure_from_humidity_ratio, 0, np.array([humidity_ratio, atmospheric_station_pressure]))
    enthalpy = np.apply_along_axis(enthalpy_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio]))
    specific_volume_moist_air = np.apply_along_axis(specific_volume_moist_air_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))
    degree_of_saturation = np.apply_along_axis(degree_of_saturation_from_humidity_ratio, 0, np.array([dry_bulb_temperature, humidity_ratio, atmospheric_station_pressure]))

    # Create DataFrame containing computed metrics
    df = pd.DataFrame(data=np.array([humidity_ratio, wet_bulb_temperature, partial_vapour_pressure_moist_air, enthalpy, specific_volume_moist_air, degree_of_saturation]).T, columns=["humidity_ratio", "wet_bulb_temperature", "partial_vapour_pressure_moist_air", "enthalpy", "specific_volume_moist_air", "degree_of_saturation"])
    print("Psychrometric calculations successful")
    return df
