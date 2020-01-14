import numpy as np

from .constants import STEFAN_BOLTZMANN_CONSTANT, KELVIN


class OpaqueMaterial(object):
    """ Material class with thermal properties

    Attributes:
        conductivity: Thermal conductivity (W/m K)
        volumetric_heat_capacity: Volumetric heat capacity (J/m^3 K)
        name: Name of the material.
    """

    def __init__(self, conductivity, volumetric_heat_capacity, specific_heat_capacity, reflectivity, emissivity, thickness=1, roughness="Medium rough"):

        # Sense check
        if emissivity <= 0 or emissivity > 1:
            raise Exception('emissivity should be > 0 and < 1')

        if roughness not in ["Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth", "Very smooth"]:
            raise Exception(
                'ground_roughness should be one of "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth", "Very smooth"')

        if reflectivity <= 0 or reflectivity > 1:
            raise Exception('reflectivity should be > 0 and < 1')

        self.conductivity = conductivity
        self.volumetric_heat_capacity = volumetric_heat_capacity # [J/kg.K] for UWG material default
        self.reflectivity = reflectivity
        self.roughness = roughness
        self.emissivity = emissivity
        self.specific_heat_capacity = specific_heat_capacity # [kJ/kg.K]
        self.thickness = thickness

    def exterior_heat_transfer_coefficient(self, local_wind_speed):
        material = {
            "Very rough": {"D": 11.58, "E": 5.89, "F": 0},
            "Rough": {"D": 12.49, "E": 4.065, "F": 0.028},
            'Medium rough': {"D": 10.79, "E": 4.192, "F": 0},
            'Medium smooth': {"D": 8.23, "E": 4, "F": -0.057},
            'Smooth': {"D": 10.22, "E": 3.1, "F": 0},
            'Very smooth': {"D": 8.23, "E": 3.33, "F": -0.036},
        }
        return material[self.roughness]["D"] + material[self.roughness][
            "E"] * local_wind_speed + material[self.roughness]["F"] * np.power(local_wind_speed, 2)

    def surface_temperature(self, dry_bulb_temperature, material_temperature, radiation_flux, local_wind_speed, kelvin=False):

        heat_transfer_coefficient = self.exterior_heat_transfer_coefficient(local_wind_speed)
        a = self.emissivity * STEFAN_BOLTZMANN_CONSTANT
        b = self.specific_heat_capacity / self.thickness + heat_transfer_coefficient
        c = -(self.specific_heat_capacity * (material_temperature + KELVIN) / self.thickness + radiation_flux * (1 - self.reflectivity) + heat_transfer_coefficient * (dry_bulb_temperature + KELVIN))

        Ts = []
        X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                  np.sqrt(3) * np.sqrt(
                              27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                                  2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
            (4 * (2 / 3) ** (1 / 3) * c) / (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                    2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
                2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
        if X.sum() > 0:
            Ts.append(X)

        X = 0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                  np.sqrt(3) * np.sqrt(
                              27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                                  2 ** (1 / 3) * 3 ** (2 / 3) * a)) + 1 / 2 * np.sqrt(-(2 * b) / (a * np.sqrt(
            (4 * (2 / 3) ** (1 / 3) * c) / (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                    2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
                2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
        if X.sum() > 0:
            Ts.append(X)

        X = -0.5 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                   np.sqrt(3) * np.sqrt(
                               27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                                   2 ** (1 / 3) * 3 ** (2 / 3) * a)) - 1 / 2 * np.sqrt((2 * b) / (a * np.sqrt(
            (4 * (2 / 3) ** (1 / 3) * c) / (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                    np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                    2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (2 ** (1 / 3) * 3 ** (
                2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (np.sqrt(3) * np.sqrt(
            27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3))
        if X.sum() > 0:
            Ts.append(X)

        X = 0.5 * np.sqrt((2 * b) / (a * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                         np.sqrt(3) * np.sqrt(
                                                     27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                         1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))) - (
                                  np.sqrt(3) * np.sqrt(
                              27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) / (
                                  2 ** (1 / 3) * 3 ** (2 / 3) * a) - (4 * (2 / 3) ** (1 / 3) * c) / (
                                  np.sqrt(3) * np.sqrt(
                              27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                  1 / 3)) - 1 / 2 * np.sqrt((4 * (2 / 3) ** (1 / 3) * c) / (
                np.sqrt(3) * np.sqrt(27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (1 / 3) + (
                                                                    np.sqrt(3) * np.sqrt(
                                                                27 * a ** 2 * b ** 4 - 256 * a ** 3 * c ** 3) + 9 * a * b ** 2) ** (
                                                                    1 / 3) / (2 ** (1 / 3) * 3 ** (2 / 3) * a))
        if X.sum() > 0:
            Ts.append(X)

        return Ts[0] if kelvin else Ts[0] - KELVIN

    def __repr__(self):
        return "Material: {0:}, k={0:}, spec vol={0:}".format(self.conductivity)
