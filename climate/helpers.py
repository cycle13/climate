import os

##########################
# Generic helper methods #
##########################


def interpret(val):
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def chunk(enumerable, n_chunks):
    for i in range(0, len(enumerable), n_chunks):
        yield enumerable[i:i + n_chunks]


def generate_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def hex_to_rgb(hex_string):
    hex_string = hex_string.lstrip('#')
    lv = len(hex_string)
    return tuple(int(hex_string[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgbList):
    if (rgbList[0] >= 1) | (rgbList[1] >= 1) | (rgbList[2] >= 1):
        return '#%02x%02x%02x' % (int(rgbList[0]), int(rgbList[1]), int(rgbList[2]))
    else:
        return '#%02x%02x%02x' % (int(rgbList[0] * 255), int(rgbList[1] * 255), int(rgbList[2] * 255))


def celsius_to_kelvin(tc):
    return tc + 273.15


#########################
# Psychrometric methods #
#########################


