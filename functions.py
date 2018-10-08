class Human:
    """Defines a single human, with physiological attributes against which comfort will be evaluated"""

    def __init__(self, age: int = 35, gender: str = "male", height: float = 1.75, weight: float = 75,
                 position: str = "standing", orientation: float = 0, clothingInsulation: float = None,
                 clothingAlbedo: float = 0.37, acclimation: bool = False, metabolicRate: float = 2.32):
        """
        :param age: Age of the human
        :param gender: Gender of the human - either "male" or "female"
        :param height: Height of the human (m)
        :param weight: Weight of the human (kg)
        :param position: Body position - either "seated" or "standing"
        :param orientation: Human facing direction - clockwise from 0 at north
        :param clothingInsulation: clo value - see Human.clothingInsulationDict()
        :param clothingAlbedo: clothing reflectivity value - see clothingAlbedoDict()
        :param acclimation: Whether the human is acclimated to evaluated conditions
        :param metabolicRate: met value - see metabolicRateDict()
        """
        self.age = age
        self.gender = gender
        self.height = height
        self.weight = weight
        self.position = position
        self.orientation = orientation
        self.clothingInsulation = clothingInsulation
        self.clothingAlbedo = clothingAlbedo
        self.acclimation = acclimation
        self.metabolicRate = metabolicRate

    def body_mass_index(self):
        """
        Return the Body Mass Index of the Human instance
        :return: Body Mass Index
        :rtype: float
        """
        return self.weight / ((self.height) ** 2)


    def body_mass_index_category(self):
        bmi = self.bodyMassIndex()
        bmi_category = None

        if self.gender == "male":
            if bmi < 17.5:
                bmi_category = "Anorexia"
            elif 17.5 <= bmi < 20.7:
                bmi_category = "Underweight"
            elif 20.7 <= bmi < 26.4:
                bmi_category = "Normal weight"
            elif 26.4 <= bmi < 27.8:
                bmi_category = "Marginally overweight"
            elif 27.8 <= bmi < 31.1:
                bmi_category = "Overweight"
            elif 31.1 <= bmi < 40:
                bmi_category = "Obese"
            elif bmi >= 40:
                bmi_category = "Extreme obesity"
        elif self.gender == "female":
            if bmi < 17.5:
                bmi_category = "Anorexia"
            elif 17.5 <= bmi < 19.1:
                bmi_category = "Underweight"
            elif 19.1 <= bmi < 25.8:
                bmi_category = "Normal weight"
            elif 25.8 <= bmi < 27.3:
                bmi_category = "Marginally overweight"
            elif 27.3 <= bmi < 32.3:
                bmi_category = "Overweight"
            elif 32.3 <= bmi < 40:
                bmi_category = "Obese"
            elif bmi >= 40:
                bmi_category = "Extreme obesity"

        return bmi_category


    def basal_metabolic_rate(self):
        if self.gender == "male":
            bmr = 9.99 * self.weight + 6.25 * (self.height * 100) - 4.92 * self.age + 5
        elif self.gender == "female":
            bmr = 9.99 * self.weight + 6.25 * (self.height * 100) - 5 * self.age - 161
        return bmr

    # Helper functions - might be temporary - who knows
    def clothing_insulation_dict(self):
        info = {
            0: "Naked!",
            0.20: "Very light summer clothes (shorts/skirt, t-shirt, slippers, no socks)",
            0.55: "Summer clothes (light trousers, short sleeves or blouse)",
            1: "Street-business suit or Typical indoor winter clothing",
            1.5: "Suit and cotton coat",
            2: "Winter suit and coat",
            2.58: "Fire-fighting clothes",
            4: "Heavy polar outfit (fur pants, coat, hood, gloves...)"
        }
        return info

    def clothing_albedo_dict(self):
        info = {
            0.21: "Dark colored (black and gray clothes)",
            0.37: "Medium colored (any clothes colors between upper two)",
            0.57: "Light colored (white and bright clothes)",
            0.95: "Protective polyethylene/aluminium suits"
        }
        return info

    def metabolic_rate_dict(self):
        info = {
            0.8: "Reclining, Sleeping",
            1: "Seated relaxed",
            1.2: "Standing at rest",
            1.4: "Car driving",
            1.5: "Graphic profession - Book Binder",
            1.6: "Standing, light activity (shopping, laboratory, light industry)",
            1.7: "Domestic work -shaving, washing and dressing",
            1.9: "Walking on the level, 2 km/h",
            2: "Standing, medium activity (shop assistant, domestic work)",
            2.2: "Building industry - Brick laying (Block of 15.3 kg)",
            2.5: "Washing dishes standing",
            2.9: "Domestic work - washing by hand, ironing, raking leaves",
            3: "Iron and steel - ramming the mold with a pneumatic hammer",
            3.1: "Building industry - forming the mold",
            3.4: "Walking on the level, 5 km/h",
            3.5: "Forestry - cutting across the grain with a one-man power saw",
            4: "Volleyball, Bicycling (15 km/h)",
            4.5: "Callisthenics",
            4.7: "Building industry - loading a wheelbarrow with stones and mortar",
            5: "Golf, Softball",
            5.5: "Gymnastics",
            6: "Aerobic Dancing, Swimming",
            6.2: "Sports - Ice skating, 18 km/h, Bicycling (20 km/h)",
            6.5: "Agriculture - digging with a spade (24 lifts/min.)",
            7: "Skiing on level (good snow, 9 km/h), Backpacking, Skating ice or roller, Basketball, Tennis",
            8: "Handball, Hockey, Racquetball, Cross County Skiing, Soccer",
            8.5: "Running 12 min/mile, Forestry - working with an axe (weight 2 kg. 33 blows/min.)",
            9.5: "Sports - Running in 15 km/h",
        }
        return info

# Method of loading weather from database
def load_weather(host, port):
    """
    Load weather from MongoDB
    :param host: Host IP address
    :type host: str
    :param port: Port number
    :type port: int
    :rtype: dict
    """
    from pymongo import MongoClient
    # TODO: Hard-code the following variables
    database = "sandp"  # Name of database
    collection = "weather"  # Name of collection
    client = MongoClient(host, port)
    db = client[database]
    col = db[collection]

    return col.find_one()

# Methods of processing climatic variables

def wind_speed_at_height(ws, h1, h2, rc=0, log=True):
    """
    Convert wind speed from one height to another using the specific method
    :param ws: wind speed (m/s)
    :type ws: float
    :param h1: original height (m)
    :type h1: float
    :param h2: target height (m)
    :type h2: float
    :param rc: roughness class corresponding to length of obstacles
    :type rc: float
    :param log: method of calculation - either log or power
    :type log: bool
    :return: temperature in Fahrenheit
    :rtype: float
    """

    import math

    roughness = {
        0: 0.0002, # Water surfaces: seas and Lakes
        0.2: 0.0005, # Inlet water
        0.5: 0.0024, # Open terrain with smooth surface, e.g. concrete, airport runways, mown grass etc.
        1: 0.03, # Open agricultural land without fences and hedges; maybe some far apart buildings and very gentle hills
        1.5: 0.055, # Agricultural land with a few buildings and 8 m high hedges separated by more than 1 km
        2: 0.1, # Agricultural land with a few buildings and 8 m high hedges separated by approx. 500 m
        2.5: 0.2, # Agricultural land with many trees, bushes and plants, or 8 m high hedges separated by approx. 250 m
        3: 0.4, # Towns, villages, agricultural land with many or high hedges, forests and very rough and uneven terrain
        3.5: 0.6, # Large towns with high buildings
        4: 1.6 # Large cities with high buildings and skyscrapers
    }

    if log:
        return ws * (math.log(h2 / roughness[rc]) / math.log(h1 / roughness[rc]))
    wind_shear_exponent = 1 / 7
    return ws * ((h2 / h1) ** wind_shear_exponent)


def ground_temperature(depth, annual_average_temperature, annual_temperature_range,
                                         days_since_coldest_day, soil_diffusivity=0.01):
    """
    Estimate ground temperature at a depth from surfaces based on air temperature
    :param depth: depth (m)
    :type depth: float
    :param annual_average_temperature: average annual air temperature (C)
    :type annual_average_temperature: float
    :param annual_temperature_range: difference between min and max annual air temperatures (C)
    :type annual_temperature_range: float
    :param days_since_coldest_day: days since the coldest day in the year
    :type days_since_coldest_day: int
    :param soil_diffusivity: soil diffusivity (cm2/s)
    :type soil_diffusivity: float

    soil_diffusivities = {
        "Rock": 0.02,
        "Wet clay": 0.015,
        "Wet sand": 0.01,
        "Dry clay": 0.002,
        "Dry sand": 0.001
    }
    """
    import math

    w = 2 * math.pi / 365
    dd = math.sqrt(2 * soil_diffusivity / w)

    return annual_average_temperature - (annual_temperature_range / 2) * math.exp(-depth / dd) * math.cos(
        (w * days_since_coldest_day) - (depth / dd))



def celsius_to_fahrenheit(Tc):
    """
    Convert temperature in Celsius to Fahrenheit
    :param Tc: temperature in Celsius
    :type Tc: float
    :return: temperature in Fahrenheit
    :rtype: float
    """
    Tf = (Tc * (9 / 5)) + 32
    return Tf


def fahrenheit_to_celsius(Tf):
    """
    Convert temperature in Fahrenheit to Celsius
    :param Tf: temperature in Fahrenheit
    :type Tf: float
    :return: temperature in Celsius
    :rtype: float
    """
    Tc = (Tf - 32) * (5 / 9)
    return Tc


def celsius_to_kelvin(Tc):
    """
    Convert temperature in Celsius to Kelvin
    :param Tc: temperature in Celsius
    :type Tc: float
    :return: temperature in Kelvin
    :rtype: float
    """
    Tk = Tc + 273.15
    return Tk


def dewpoint_temperature(Ta, rh):
    # TODO: Add check that rh is between 0.0 and 100.0 for validity - if not this may not be as accurate!
    """
    Calculate dew-point temperature from air temperature and relative humidity using MET4 and MET4A method
    for Tc in range:  0C < Tc < 60C and Tdp in range:  0C < Tdp < 50C, uncertainty is +/- 0.4C
    :param Ta: temperature (Celsius)
    :type Ta: float
    :param rh: relative humidity (%)
    :type rh: float
    :return: dew-point temperature (Celsius)
    :rtype: float
    """

    import math
    import warnings

    if not 0 <= rh <= 100:
        warnings.warn("The value input for Relative Humidity is outside the range suitable for accurate Dewpoint estimation")

    a = 17.27
    b = 237.7
    rh = rh / 100
    Tdp = (b * (((a * Ta) / (b + Ta)) + math.log(rh))) / (a - ((( a * Ta) / (b + Ta)) + math.log(rh)))

    return Tdp


def humidity_index(Ta, Tdp):
    """
    Calulate the humidity index (humidex) from air temperature and relative humidity -  see https://web.archive.org/web/20130627223738/http://climate.weatheroffice.gc.ca/prods_servs/normals_documentation_e.html
    :param Ta: temperature (Celcius)
    :type Ta: float
    :param Tdp: dew-point temperature (Celsius)
    :type Tdp: float
    :return: humidity index (%), humidity index category (0-4), comfortable
    :rtype: Tuple[float, int, bool]
    """

    import math

    dewpoint_k = celsius_to_kelvin(Tdp)

    e = 6.11 * math.exp(5417.7530 * ((1 / 273.16) - (1 / dewpoint_k)))
    h = (0.5555) * (e - 10.0)
    humidity_index = Ta + h

    humidity_index_categories = {
        0: "Little to no discomfort",
        1: "Some discomfort",
        2: "Great discomfort; avoid exertion",
        3: "Dangerous; heat stroke quite possible"
    }

    humidity_index_category = None
    comfortable = None

    if humidity_index < 30:
        # Little to no discomfort
        humidity_index_category = 0
        comfortable = True
    elif 30 <= humidity_index < 39:
        # Some discomfort
        humidity_index_category = 1
        comfortable = False
    elif 40 <= humidity_index < 45:
        # Great discomfort; avoid exertion
        humidity_index_category = 2
        comfortable = False
    elif humidity_index >= 45:
        # Dangerous; heat stroke quite possible
        humidity_index_category = 3
        comfortable = False

    return humidity_index, humidity_index_category, comfortable


def heat_index(Ta, rh):
    """
    Compute the Heat Index in accordance with the regression defined in the National Weather Service (NWS) Technical
    Attachment (SR 90-23) - for details see https://www.weather.gov/safety/heat-index
    :param Ta: temperature in Celcius
    :type Ta: float
    :param rh: relative humidity
    :type rh: float
    :return: heat index in Celcius, heat index category (0-4), comfortable
    :rtype: Tuple[float, int, bool]
    """

    import math

    Tf = celsius_to_fahrenheit(Ta)

    if Tf < 80:
        heat_index_fahrenheit = 0.5 * (Tf + 61.0 + ((Tf - 68.0) * 1.2) + (rh * 0.094))
    else:
        heat_index_fahrenheit = -42.379 + 2.04901523 * Tf + 10.14333127 * rh - 0.22475541 * Tf * rh - 6.83783 * (
                    10 ** (-3)) * (Tf ** (2)) - 5.481717 * (10 ** (-2)) * (rh ** (2)) + 1.22874 * (10 ** (-3)) * (
                                            Tf ** (2)) * (rh) + 8.5282 * (10 ** (-4)) * (Tf) * (rh ** (2)) - 1.99 * (
                                            10 ** (-6)) * (Tf ** (2)) * (rh ** (2))
        if (Tf >= 80) and (Tf <= 112) and (rh < 13):
            adjust = ((13 - rh) / 4) * math.sqrt((17 - abs(Tf - 95)) / 17)
            heat_index_fahrenheit = heat_index_fahrenheit - adjust
        elif (Tf >= 80) and (Tf <= 87) and (rh > 85):
            adjust = ((rh - 85) / 10) * ((87 - Tf) / 5)
            heat_index_fahrenheit = heat_index_fahrenheit + adjust

    heat_index_categories = {
        0: None,
        1: "Caution",
        2: "Extreme caution",
        3: "Danger",
        4: "Extreme danger"
    }
    heat_index_category = None
    comfortable = None

    if heat_index_fahrenheit < 80:
        heat_index_category = 0
        comfortable = True
    elif 80 <= heat_index_fahrenheit < 90:
        heat_index_category = 1
        comfortable = False
    elif 90 <= heat_index_fahrenheit < 105:
        heat_index_category = 2
        comfortable = False
    elif 105 <= heat_index_fahrenheit < 130:
        heat_index_category = 3
        comfortable = False
    elif heat_index_fahrenheit >= 130:
        heat_index_category = 4
        comfortable = False

    heat_index = fahrenheit_to_celsius(heat_index_fahrenheit)

    return heat_index, heat_index_category, comfortable


def discomfort_index(Ta, rh):
    """
    Compute the Heat Index (Thom's Index) in accordance with the regression defined in Thom, E.C. (1959): The discomfort index. Weather wise, 12: 5760
    :param Ta: temperature in Celcius
    :type Ta: float
    :param rh: relative humidity
    :type rh: float
    :return: discomfort index in Celcius, discomfort index category (0-9), comfortable
    :rtype: Tuple[float, int, bool]
    """

    discomfort_index = Ta - (0.55 - 0.0055 * rh) * (Ta - 14.5)

    discomfort_index_categories = {
        0: "Hyperglacial",
        1: "Glacial",
        2: "Extremely cold",
        3: "Very cold",
        4: "Cold",
        5: "Cool",
        6: "Comfortable",
        7: "Hot",
        8: "Very hot",
        9: "Torrid"
    }

    # Categories from Kyle WJ (1994) The human bioclimate of Hong-Kong.
    discomfort_index_category = None
    comfortable = None
    if discomfort_index < -40:
        discomfort_index_category = 0
        comfortable = False
    elif -40 <= discomfort_index < -20:
        discomfort_index_category = 1
        comfortable = False
    elif -20 <= discomfort_index < -10:
        discomfort_index_category = 2
        comfortable = False
    elif -10 <= discomfort_index < -1.8:
        discomfort_index_category = 3
        comfortable = False
    elif -1.8 <= discomfort_index < 13:
        discomfort_index_category = 4
        comfortable = False
    elif 13 <= discomfort_index < 15:
        discomfort_index_category = 5
        comfortable = False
    elif 15 <= discomfort_index < 20:
        discomfort_index_category = 6
        comfortable = True
    elif 20 <= discomfort_index < 26.5:
        discomfort_index_category = 7
        comfortable = False
    elif 26.5 <= discomfort_index < 30:
        discomfort_index_category = 8
        comfortable = False
    elif discomfort_index >= 30:
        discomfort_index_category = 9
        comfortable = False

    return discomfort_index, discomfort_index_category, comfortable


def wind_chill_temperature(Ta, ws):
    """
    Compute the Wind Chill Index, based on the Joint Action Group for Temperature Indices (JAG/TI) methodology
    :param Ta: temperature in Celcius
    :type Ta: float
    :param ws: wind speed (m/s)
    :type ws: float
    :return: wind chill index, wind chill index category (0-6), comfortable
    :rtype: Tuple[float, int, bool]
    """

    ws_km_h = ws * 3.6   # convert m/s to km/h wind speed

    wind_chill_temperature = 13.12 + 0.6215 * Ta - 11.37 * (ws_km_h ** 0.16) + 0.3965 * Ta * (ws_km_h ** 0.16)

    wind_chill_index_categories = {
        0: "No risk",
        1: "Low risk: Slight increase in discomfort",
        2: "Moderate risk: Risk of hypothermia and frostbite",
        3: "High risk: High risk of frostbite and hypothermia. Exposed skin can freeze in 10 to 30 minutes",
        4: "Very high risk: Very high risk of frostbite and hypothermia. Exposed skin can freeze in 5 to 10 minutes",
        5: "Severe risk: Severe risk of frostbite and hypothermia. Exposed skin can freeze in 2 to 5 minutes",
        6: "Extreme risk: Outdoor conditions are hazardous. Exposed skin can freeze in less than 2 minutes"
    }

    wind_chill_index_category = None
    comfortable = None

    if wind_chill_temperature >= 0:
        wind_chill_index_category = 0
        comfortable = True
    elif 0 > wind_chill_temperature >= -9:
        wind_chill_index_category = 1
        comfortable = False
    elif -9 > wind_chill_temperature >= -27:
        wind_chill_index_category = 2
        comfortable = False
    elif -27 > wind_chill_temperature >= -39:
        wind_chill_index_category = 3
        comfortable = False
    elif -39 > wind_chill_temperature >= -47:
        wind_chill_index_category = 4
        comfortable = False
    elif -47 > wind_chill_temperature >= -54:
        wind_chill_index_category = 5
        comfortable = False
    elif wind_chill_temperature < -54:
        wind_chill_index_category = 6
        comfortable = False

    return wind_chill_temperature, wind_chill_index_category, comfortable


def wet_bulb_globe_temperature_outdoors(Ta, Tr, vp, ws):
    """
    Compute the wet-bulb globe temperature from air temperature, mean radiant temperature and vapour pressure
    :param Ta: temperature in Celcius
    :type Ta: float
    :param Tr: mean radiant temperature (Celsius)
    :type Tr: float
    :param vp: vapour pressure (Pa)
    :type vp: float
    :param ws: wind speed (m/s)
    :type ws: float
    :return: WBGT, WBGT category (0-4), comfortable
    :rtype: Tuple[float, int, bool]
    """

    wet_bulb_temperature = 1.885 + 0.3704 * Ta + 0.4492 * (vp / 1000)
    globe_temperature = 2.098 - 2.561 * ws + 0.5957 * Ta + 0.4017 * Tr
    wet_bulb_globe_temperature = 0.7 * wet_bulb_temperature + 0.2 * globe_temperature + 0.1 * Ta

    wet_bulb_globe_temperature_categories = {
        0: "No heat stress",
        1: "Minor heat stress",
        2: "Heat stress",
        3: "High risk: High risk of frostbite and hypothermia. Exposed skin can freeze in 10 to 30 minutes",
        4: "Very high risk: Very high risk of frostbite and hypothermia. Exposed skin can freeze in 5 to 10 minutes",
        5: "Severe risk: Severe risk of frostbite and hypothermia. Exposed skin can freeze in 2 to 5 minutes"
    }

    wet_bulb_globe_temperature_category = None
    comfortable = None

    if wet_bulb_globe_temperature < 18:
        wet_bulb_globe_temperature_category = 0
        comfortable = True
    elif 18 <= wet_bulb_globe_temperature < 23:
        wet_bulb_globe_temperature_category = 1
        comfortable = False
    elif 23 <= wet_bulb_globe_temperature < 28:
        wet_bulb_globe_temperature_category = 2
        comfortable = False
    elif 28 <= wet_bulb_globe_temperature < 30:
        wet_bulb_globe_temperature_category = 3
        comfortable = False
    elif wet_bulb_globe_temperature >= 30:
        wet_bulb_globe_temperature_category = 4
        comfortable = False

    return wet_bulb_globe_temperature, wet_bulb_globe_temperature_category, comfortable


def effective_temperature(Ta, rh, ws, SR, ac):
    """
    Compute the effective temperature
    :param Ta: temperature (C)
    :type Ta: float
    :param ac: clothing insulation (%)
    :type ac: float
    :param SR: solar radiation (Wh/m2)
    :type SR: float
    :param rh: relative humidity (%)
    :type rh: float
    :param ws: wind speed (m/s)
    :type ws: float
    :return: effective temperature, effective temperature category (0-6), comfortable
    :rtype: Tuple[float, int, bool]
    """

    effective_temperature = None

    if ws <= 0.2:
        effective_temperature = Ta - 0.4 * (Ta - 10) * (1 - rh / 100)  # formula by Missenard
    elif ws > 0.2:
        effective_temperature = 37 - ((37 - Ta) / (0.68 - (0.0014 * rh) + (1 / (1.76 + 1.4 * (ws ** 0.75))))) - (
                    0.29 * Ta * (1 - 0.01 * rh))  # modified formula by Gregorczuk (WMO, 1972; Hentschel, 1987)

    # Radiative-effective temperature
    radiative_effective_temperature = effective_temperature + ((1 - 0.01 * ac) * SR) * (
                (0.0155 - 0.00025 * effective_temperature) - (0.0043 - 0.00011 * effective_temperature))

    effective_temperature_category = None
    comfortable = None

    if radiative_effective_temperature < 1:
        effective_temperature_category = 0
        comfortable = False
    elif radiative_effective_temperature >= 1 and effective_temperature < 9:
        effective_temperature_category = 1
        comfortable = False
    elif radiative_effective_temperature >= 9 and effective_temperature < 17:
        effective_temperature_category = 2
        comfortable = False
    elif radiative_effective_temperature >= 17 and effective_temperature < 21:
        effective_temperature_category = 3
        comfortable = False
    elif radiative_effective_temperature >= 21 and effective_temperature < 23:
        effective_temperature_category = 4
        comfortable = True
    elif radiative_effective_temperature >= 23 and effective_temperature < 27:
        effective_temperature_category = 5
        comfortable = False
    elif radiative_effective_temperature >= 27:
        effective_temperature_category = 6
        comfortable = False

    return radiative_effective_temperature, effective_temperature_category, comfortable


# def comfPierceSET(Ta, Tr, ws, rh, met, clo, wme):
#     """
#     Compute the standard effective temperature
#     :param Ta: temperature (C)
#     :type Ta: float
#     :param Tr: mean radiant temperature (C)
#     :type Tr: float
#     :param ws: wind speed (m/s)
#     :type ws: float
#     :param rh: relative humidity (%)
#     :type rh: float
#     :param met: metabolic rate (met)
#     :type met: float
#     :param clo: clo value ()
#     :type clo: float
#     :param wme: metabolic rate (met)
#     :type wme: float
#     :return: standard effective temperature
#     :rtype: float
#     """
#
#     metabolic_rates = {
#         "Sleeping": 0.7,
#         "Reclining": 0.8,
#         "Sitting": 1.0,
#         "Typing": 1.1,
#         "Standing": 1.2,
#         "Driving": 1.5,
#         "Cooking": 1.8,
#         "Walking": 1.7,
#         "Walking 2mph": 2.0,
#         "Lifting 10lbs": 2.1,
#         "Walking 3mph": 2.6,
#         "House Cleaning": 2.7,
#         "Basketball": 3,
#         "Dancing": 3.4,
#         "Walking 4mph": 3.8,
#         "Lifting 100lbs": 4.0,
#         "Shoveling": 4.4,
#         "Running 9mph": 9.5
# }
#
#     #Function to find the saturation vapour pressure, used frequently throughout the comfPierceSET function.
#     def findSaturatedVaporPressureTorr(T):
#         #calculates Saturated Vapor Pressure (Torr) at Temperature T  (C)
#         return np.exp(18.6686 - 4030.183 / (T + 235.0))
#
#     #Key initial variables.
#     VaporPressure = (rh * findSaturatedVaporPressureTorr(Ta)) / 100
#     AirVelocity = max(ws, 0.1)
#     KCLO = 0.25
#     BODYWEIGHT = 69.9
#     BODYSURFACEAREA = 1.8258
#     METFACTOR = 58.2
#     SBC = 0.000000056697 # Stefan-Boltzmann constant (W/m2K4)
#     CSW = 170
#     CDIL = 120
#     CSTR = 0.5
#
#     TempSkinNeutral = 33.7 #setpoint (neutral) value for Tsk
#     TempCoreNeutral = 36.8 #setpoint value for Tcr
#     TempBodyNeutral = 36.49 #setpoint for Tb (.1*TempSkinNeutral + .9*TempCoreNeutral)
#     SkinBloodFlowNeutral = 6.3 #neutral value for SkinBloodFlow
#
#     #INITIAL VALUES - start of 1st experiment
#     TempSkin = TempSkinNeutral
#     TempCore = TempCoreNeutral
#     SkinBloodFlow = SkinBloodFlowNeutral
#     MSHIV = 0.0
#     ALFA = 0.1
#     ESK = 0.1 * met
#
#     #Start new experiment here (for graded experiments)
#     #UNIT CONVERSIONS (from input variables)
#
#     p = 101325.0 / 1000 # This variable is the pressure of the atmosphere in kPa and was taken from the psychrometrics.js file of the CBE comfort tool.
#
#     PressureInAtmospheres = p * 0.009869
#     LTIME = 60
#     TIMEH = LTIME / 60.0
#     RCL = 0.155 * clo
#     #AdjustICL(RCL, Conditions);  TH: I don't think this is used in the software
#
#     FACL = 1.0 + 0.15 * clo #% INCREASE IN BODY SURFACE AREA DUE TO CLOTHING
#     LR = 2.2 / PressureInAtmospheres #Lewis Relation is 2.2 at sea level
#     RM = met * METFACTOR
#     M = met * METFACTOR
#
#     if clo <= 0:
#         WCRIT = 0.38 * pow(AirVelocity, -0.29)
#         ICL = 1.0
#     else:
#         WCRIT = 0.59 * pow(AirVelocity, -0.08)
#         ICL = 0.45
#
#     CHC = 3.0 * pow(PressureInAtmospheres, 0.53)
#     CHCV = 8.600001 * pow((AirVelocity * PressureInAtmospheres), 0.53)
#     CHC = max(CHC, CHCV)
#
#     #initial estimate of Tcl
#     CHR = 4.7
#     CTC = CHR + CHC
#     RA = 1.0 / (FACL * CTC) #resistance of air layer to dry heat transfer
#     TOP = (CHR * Tr + CHC * Ta) / CTC
#     TCL = TOP + (TempSkin - TOP) / (CTC * (RA + RCL))
#
#     # ========================  BEGIN ITERATION
#     #
#     # Tcl and CHR are solved iteratively using: H(Tsk - To) = CTC(Tcl - To),
#     # where H = 1/(Ra + Rcl) and Ra = 1/Facl*CTC
#
#     TCL_OLD = TCL
#     TIME = range(LTIME)
#     flag = True
#     for TIM in TIME:
#         if flag == True:
#             while abs(TCL - TCL_OLD) > 0.01:
#                 TCL_OLD = TCL
#                 CHR = 4.0 * SBC * pow(((TCL + Tr) / 2.0 + 273.15), 3.0) * 0.72
#                 CTC = CHR + CHC
#                 RA = 1.0 / (FACL * CTC) #resistance of air layer to dry heat transfer
#                 TOP = (CHR * Tr + CHC * Ta) / CTC
#                 TCL = (RA * TempSkin + RCL * TOP) / (RA + RCL)
#         flag = False
#         DRY = (TempSkin - TOP) / (RA + RCL)
#         HFCS = (TempCore - TempSkin) * (5.28 + 1.163 * SkinBloodFlow)
#         ERES = 0.0023 * M * (44.0 - VaporPressure)
#         CRES = 0.0014 * M * (34.0 - Ta)
#         SCR = M - HFCS - ERES - CRES - wme
#         SSK = HFCS - DRY - ESK
#         TCSK = 0.97 * ALFA * BODYWEIGHT
#         TCCR = 0.97 * (1 - ALFA) * BODYWEIGHT
#         DTSK = (SSK * BODYSURFACEAREA) / (TCSK * 60.0)# //deg C per minute
#         DTCR = SCR * BODYSURFACEAREA / (TCCR * 60.0)# //deg C per minute
#         TempSkin = TempSkin + DTSK
#         TempCore = TempCore + DTCR
#         TB = ALFA * TempSkin + (1 - ALFA) * TempCore
#         SKSIG = TempSkin - TempSkinNeutral
#         WARMS = (SKSIG > 0) * SKSIG
#         COLDS = ((-1.0 * SKSIG) > 0) * (-1.0 * SKSIG)
#         CRSIG = (TempCore - TempCoreNeutral)
#         WARMC = (CRSIG > 0) * CRSIG
#         COLDC = ((-1.0 * CRSIG) > 0) * (-1.0 * CRSIG)
#         BDSIG = TB - TempBodyNeutral
#         WARMB = (BDSIG > 0) * BDSIG
#         COLDB = ((-1.0 * BDSIG) > 0) * (-1.0 * BDSIG)
#         SkinBloodFlow = (SkinBloodFlowNeutral + CDIL * WARMC) / (1 + CSTR * COLDS)
#         if SkinBloodFlow > 90.0: SkinBloodFlow = 90.0
#         if SkinBloodFlow < 0.5: SkinBloodFlow = 0.5
#         REGSW = CSW * WARMB * np.exp(WARMS / 10.7)
#         if REGSW > 500.0: REGSW = 500.0
#         ERSW = 0.68 * REGSW
#         REA = 1.0 / (LR * FACL * CHC) #evaporative resistance of air layer
#         RECL = RCL / (LR * ICL) #evaporative resistance of clothing (icl=.45)
#         EMAX = (findSaturatedVaporPressureTorr(TempSkin) - VaporPressure) / (REA + RECL)
#         PRSW = ERSW / EMAX
#         PWET = 0.06 + 0.94 * PRSW
#         EDIF = PWET * EMAX - ERSW
#         ESK = ERSW + EDIF
#         if PWET > WCRIT:
#             PWET = WCRIT
#             PRSW = WCRIT / 0.94
#             ERSW = PRSW * EMAX
#             EDIF = 0.06 * (1.0 - PRSW) * EMAX
#             ESK = ERSW + EDIF
#         if EMAX < 0:
#             EDIF = 0
#             ERSW = 0
#             PWET = WCRIT
#             PRSW = WCRIT
#             ESK = EMAX
#         ESK = ERSW + EDIF
#         MSHIV = 19.4 * COLDS * COLDC
#         M = RM + MSHIV
#         ALFA = 0.0417737 + 0.7451833 / (SkinBloodFlow + .585417)
#
#     #Define new heat flow terms, coeffs, and abbreviations
#     STORE = M - wme - CRES - ERES - DRY - ESK #rate of body heat storage
#     HSK = DRY + ESK #total heat loss from skin
#     RN = M - wme #net metabolic heat production
#     ECOMF = 0.42 * (RN - (1 * METFACTOR))
#     if ECOMF < 0.0: ECOMF = 0.0 #from Fanger
#     EREQ = RN - ERES - CRES - DRY
#     EMAX = EMAX * WCRIT
#     HD = 1.0 / (RA + RCL)
#     HE = 1.0 / (REA + RECL)
#     W = PWET
#     PSSK = findSaturatedVaporPressureTorr(TempSkin)
#     #Definition of ASHRAE standard environment... denoted "S"
#     CHRS = CHR
#     if met < 0.85:
#         CHCS = 3.0
#     else:
#         CHCS = 5.66 * pow((met - 0.85), 0.39)
#         if CHCS < 3.0: CHCS = 3.0
#
#     CTCS = CHCS + CHRS
#     RCLOS = 1.52 / ((met - wme / METFACTOR) + 0.6944) - 0.1835
#     RCLS = 0.155 * RCLOS
#     FACLS = 1.0 + KCLO * RCLOS
#     FCLS = 1.0 / (1.0 + 0.155 * FACLS * CTCS * RCLOS)
#     IMS = 0.45
#     ICLS = IMS * CHCS / CTCS * (1 - FCLS) / (CHCS / CTCS - FCLS * IMS)
#     RAS = 1.0 / (FACLS * CTCS)
#     REAS = 1.0 / (LR * FACLS * CHCS)
#     RECLS = RCLS / (LR * ICLS)
#     HD_S = 1.0 / (RAS + RCLS)
#     HE_S = 1.0 / (REAS + RECLS)
#
#     #SET* (standardized humidity, clo, Pb, and CHC) determined using Newton's iterative solution
#     #FNERRS is defined in the GENERAL SETUP section above
#
#     DELTA = .0001
#     dx = 100.0
#     X_OLD = TempSkin - HSK / HD_S #lower bound for SET
#     while abs(dx) > .01:
#         ERR1 = (HSK - HD_S * (TempSkin - X_OLD) - W * HE_S * (PSSK - 0.5 * findSaturatedVaporPressureTorr(X_OLD)))
#         ERR2 = (HSK - HD_S * (TempSkin - (X_OLD + DELTA)) - W * HE_S * (PSSK - 0.5 * findSaturatedVaporPressureTorr((X_OLD + DELTA))))
#         X = X_OLD - DELTA * ERR1 / (ERR2 - ERR1)
#         dx = X - X_OLD
#         X_OLD = X
#
#     return X