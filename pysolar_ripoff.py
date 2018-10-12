from datetime import utctimetuple

def tm_yday(d):
    return d.utctimetuple().tm_yday

def get_azimuth_fast(latitude_deg, longitude_deg, when):
    # expect 230 degrees for solar.get_azimuth(42.364908,-71.112828,datetime.datetime(2007, 2, 18, 20, 18, 0, 0))
    day = tm_yday(when)
    declination_rad = math.radians(get_declination(day))
    latitude_rad = math.radians(latitude_deg)
    hour_angle_rad = math.radians(get_hour_angle(when, longitude_deg))
    altitude_rad = math.radians(get_altitude_fast(latitude_deg, longitude_deg, when))

    azimuth_rad = math.asin(-math.cos(declination_rad) * math.sin(hour_angle_rad) / math.cos(altitude_rad))

    return math.where(math.cos(hour_angle_rad) * math.tan(latitude_rad) >= math.tan(declination_rad),
                      (180 - math.degrees(azimuth_rad)),
                      math.degrees(azimuth_rad) + 360 * (azimuth_rad < 0)
                      )