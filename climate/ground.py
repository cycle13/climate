import pandas as pd
import numpy as np

from .helpers import chunk




# def ground_temperature_at_depth(depth, annual_average_temperature, annual_temperature_range, days_since_coldest_day, soil_diffusivity="Dry clay"):
#     soil_diffusivities = {
#         "Rock": 0.02,
#         "Wet clay": 0.015,
#         "Wet sand": 0.01,
#         "Dry clay": 0.002,
#         "Dry sand": 0.001
#     }
#
#     w = 2 * np.pi / 365
#     dd = np.sqrt(2 * soil_diffusivities[soil_diffusivity] / w)
#
#     return annual_average_temperature - (annual_temperature_range / 2) * np.exp(-depth / dd) * np.cos((w * days_since_coldest_day) - (depth / dd))
#
#
# def annual_ground_temperature_at_depth(self, depth, soil_diffusivity="Dry clay"):
#     annual_average_temperature = self.dry_bulb_temperature.mean()
#     annual_temperature_range = self.dry_bulb_temperature.max() - self.dry_bulb_temperature.min()
#     days_since_coldest_day = np.array([i if i > 0 else i + 365 for i in (self.index - self.dry_bulb_temperature.resample("1D").mean().idxmin()).total_seconds() / 86400])
#     self.ground_temperature_calculated = pd.Series(index=self.index, name="ground_temperature_{0:0.0f}_calculated".format(depth * 1000), data=ground_temperature_at_depth(depth, annual_average_temperature, annual_temperature_range, days_since_coldest_day, soil_diffusivity))
#     print("Ground temperature approximation successful")
#     return self.ground_temperature_calculated
#
