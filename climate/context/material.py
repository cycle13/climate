import numpy as np
import uuid

from ..constants import STEFAN_BOLTZMANN_CONSTANT, KELVIN

class Material(object):
    def __init__(self, thickness: float = 0.1, roughness: str = "MediumRough", conductivity: float = 0.75,
                 density: float = 2243, specific_heat_capacity: float = 900, reflectivity: float = 0.05, transmissivity: float = 0,
                 emissivity: float = 0.93):

        self.guid = str(uuid.uuid4())
        self.thickness = thickness
        self.roughness = roughness
        self.conductivity = conductivity
        self.density = density
        self.specific_heat_capacity = specific_heat_capacity
        self.reflectivity = reflectivity
        self.transmissivity = transmissivity
        self.absorptivity = 1 - reflectivity - transmissivity
        self.emissivity = emissivity  # in E+ as thermal absorptance

        self.sense_check()

    def sense_check(self):
        # Sense check
        roughnesses = ["VeryRough", "Rough", "MediumRough", "MediumSmooth", "Smooth", "VerySmooth"]
        if self.roughness not in roughnesses:
            raise Exception('Roughness should be one of {}'.format(str(roughnesses)))

        if not 0.003 < self.thickness <= 1:
            raise Exception('Thickness should be > 0.003 and <= 1')

        if not 0.01 < self.conductivity <= 5:
            raise Exception('Conductivity should be > 0.01 and <= 5')

        if not 10 < self.density <= 8000:
            raise Exception('Density should be > 10 and <= 8000')

        if not 200 < self.specific_heat_capacity <= 2000:
            raise Exception('Density should be > 200 and <= 2000')

        if not 0 < self.reflectivity < 1:
            raise Exception('Reflectivity should be > 0 and < 1')

        if not 0 <= self.transmissivity < 1:
            raise Exception('Transmissivity should be > 0 and < 1')

        if not 0 < self.emissivity < 1:
            raise Exception('Emissivity should be > 0 and < 1')

    def exterior_heat_transfer_coefficient(self, local_wind_speed):
        return exterior_heat_transfer_coefficient(self.roughness, local_wind_speed)

    def surface_temperature(self, dry_bulb_temperature, local_wind_speed, radiation_flux, material_temperature,  kelvin=False):
        return surface_temperature(dry_bulb_temperature, local_wind_speed, radiation_flux, material_temperature, self.reflectivity, self.emissivity, self.specific_heat_capacity, self.thickness, self.roughness, kelvin=kelvin)

    def to_rad_string(self):
        return rad_string_plastic(id=self.guid, reflectivity=self.reflectivity) if self.transmissivity == 0 else rad_string_trans(id=self.guid, transmissivity=self.transmissivity, specularity=1.52)

    def __repr__(self):
        return_string = "Material: \n"
        for k, v in self.__dict__.items():
            return_string += "- {}: {}\n".format(k.title(), v)
        return return_string


def rad_string_plastic(id: str=str(uuid.uuid4()), reflectivity: float=0.5):
    return "void plastic {0:}\n0\n0\n5 {1:} {1:} {1:} 0.0 0.0".format(id, reflectivity)


def rad_string_trans(id: str=str(uuid.uuid4()), transmissivity: float=0.5, specularity: float=1.52):
    return "void glass {0:}\n0\n0\n4 {1:} {1:} {1:} {2:}".format(id, transmissivity, specularity)


def exterior_heat_transfer_coefficient(material_roughness, local_wind_speed):
    """
    The simple algorithm uses surface roughness and local surface windspeed to calculate the exterior heat transfer coefficient.

    Note that the simple correlation yields a combined convection and radiation heat transfer coefficient.

    Material roughness correlations are taken from Figure 1, Page 22.4, ASHRAE Handbook of Fundamentals (ASHRAE 1989).

    Source: http://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-020.html for further details.

    Parameters
    ----------
    material_roughness : string
        A surface roughness index. Acceptable values are: "Very rough", "Rough", "Medium rough", "Medium smooth", "Smooth" and "Very smooth"
    local_wind_speed : float
        Local wind speed calculated at the height above ground of the surface centroid

    Returns
    -------
    heat_transfer_coefficient
        The surface heat transfer coefficient for the variables passed
    """
    material = {
        "VeryRough": {"D": 11.58, "E": 5.89, "F": 0},
        "Rough": {"D": 12.49, "E": 4.065, "F": 0.028},
        'MediumRough': {"D": 10.79, "E": 4.192, "F": 0},
        'MediumSmooth': {"D": 8.23, "E": 4, "F": -0.057},
        'Smooth': {"D": 10.22, "E": 3.1, "F": 0},
        'VerySmooth': {"D": 8.23, "E": 3.33, "F": -0.036},
    }
    return material[material_roughness]["D"] + material[material_roughness][
        "E"] * local_wind_speed + material[material_roughness]["F"] * np.power(local_wind_speed, 2)


def surface_temperature(dry_bulb_temperature, local_wind_speed, radiation_flux, material_temperature, material_reflectivity, material_emissivity, material_specific_heat_capacity, material_thickness, material_roughness, kelvin=False):

    heat_transfer_coefficient = exterior_heat_transfer_coefficient(material_roughness, local_wind_speed)
    a = material_emissivity * STEFAN_BOLTZMANN_CONSTANT
    b = material_specific_heat_capacity / material_thickness + heat_transfer_coefficient
    c = -(material_specific_heat_capacity * (material_temperature + KELVIN) / material_thickness + radiation_flux * (1 - material_reflectivity) + heat_transfer_coefficient * (dry_bulb_temperature + KELVIN))

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