"""Material class"""


class Material(object):
    """ uwg Material

    Attributes:
        thermalCond: Thermal conductivity (W/m K)
        volHeat: Volumetric heat capacity (J/m^3 K)
        name: Name of the material.
    """

    def __init__(self, thermalCond, volHeat, name='noname'):
        self._name = name  # purely for internal purpose
        self.thermalCond = thermalCond
        self.volHeat = volHeat

    def __repr__(self):
        return "Material: {0:}, k={1:}, spec vol={2:}".format(self._name, self.thermalCond, self.volHeat)
