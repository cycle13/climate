import numpy as np
import pandas as pd

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


def pedestrian_wind_speed(self, source_wind_height=10, target_wind_height=1.5, terrain_roughness="Airport runway areas",
                          log_method=True):
    self.pedestrian_wind_speed = pd.Series(name="pedestrian_wind_speed", index=self.index,
                                           data=wind_speed_at_height(self.wind_speed, source_wind_height,
                                                                     target_wind_height, terrain_roughness, log_method))
    print("Wind speed translation successful")
    return self.pedestrian_wind_speed
