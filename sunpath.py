"""
This code was taken from Ladybug-Tools, https://github.com/ladybug-tools/ladybug

Currently, it's messy, contains a lot of superfluous functionality and needs simplifying and sorting into something
much more lightweight and function specific (i.e. calculating the solar altitude and azimuth only).

A lot of this code should be replaced with more accurate and detailed code from PySolar!

Use it by importing the file, then running ...

import sunpath
s = sunpath.Sunpath(latitude=latitude, longitude=longitude, time_zone=time_zone, north_angle=0, daylight_saving_period=False)
dt = sunpath.DateTime(month=1, day=1, hour=1, minute=0)
s.calculate_sun_from_date_time(dt).azimuth
s.calculate_sun_from_date_time(dt).altitude

"""


from datetime import datetime
import re
import math
from collections import namedtuple
import sys
if (sys.version_info > (3, 0)):
    from euclid3 import Vector3
    xrange = range
else:
    from euclid import Vector3

class DateTime(datetime):
    """Create Ladybug Date time.
    Args:
        month: A value for month between 1-12 (Defualt: 1).
        day: A value for day between 1-31 (Defualt: 1).
        hour: A value for hour between 0-23 (Defualt: 0).
        minute: A value for month between 0-59 (Defualt: 0).
        leap_year: A boolean to indicate if datetime is for a leap year
            (Default: False).
    """

    __slots__ = ()

    def __new__(cls, month=1, day=1, hour=0, minute=0, leap_year=False):
        """Create Ladybug datetime.
        Args:
            month: A value for month between 1-12 (Defualt: 1).
            day: A value for day between 1-31 (Defualt: 1).
            hour: A value for hour between 0-23 (Defualt: 0).
            minute: A value for month between 0-59 (Defualt: 0).
            leap_year: A boolean to indicate if datetime is for a leap year
                (Default: False).
        """
        year = 2016 if leap_year else 2017
        hour, minute = cls._calculate_hour_and_minute(hour + minute / 60.0)
        try:
            return datetime.__new__(cls, year, month, day, hour, minute)
        except ValueError as e:
            raise ValueError("{}:\n\t({}/{}@{}:{})(m/d@h:m)".format(
                e, month, day, hour, minute
            ))

    @classmethod
    def from_json(cls, data):
        """Creat datetime from a dictionary.
        Args:
            data: {
                'month': A value for month between 1-12. (Defualt: 1)
                'day': A value for day between 1-31. (Defualt: 1)
                'hour': A value for hour between 0-23. (Defualt: 0)
                'minute': A value for month between 0-59. (Defualt: 0)
            }
        """
        if 'month' not in data:
            data['month'] = 1

        if 'day' not in data:
            data['day'] = 1

        if 'hour' not in data:
            data['hour'] = 0

        if 'minute' not in data:
            data['minute'] = 0

        if 'year' not in data:
            data['year'] = 2017

        leap_year = True if int(data['year']) == 2016 else False
        return cls(data['month'], data['day'], data['hour'], data['minute'], leap_year)

    @classmethod
    def from_hoy(cls, hoy, leap_year=False):
        """Create Ladybug Datetime from an hour of the year.
        Args:
            hoy: A float value 0 <= and < 8760
        """
        return cls.from_moy(round(hoy * 60), leap_year)

    @classmethod
    def from_moy(cls, moy, leap_year=False):
        """Create Ladybug Datetime from a minute of the year.
        Args:
            moy: An integer value 0 <= and < 525600
        """
        if not leap_year:
            num_of_minutes_until_month = (0, 44640, 84960, 129600, 172800, 217440,
                                          260640, 305280, 349920, 393120, 437760,
                                          480960, 525600)
        else:
            num_of_minutes_until_month = (0, 44640, 84960 + 1440, 129600 + 1440,
                                          172800 + 1440, 217440 + 1440, 260640 + 1440,
                                          305280 + 1440, 349920 + 1440, 393120 + 1440,
                                          437760 + 1440, 480960 + 1440, 525600 + 1440)
        # find month
        for monthCount in range(12):
            if int(moy) < num_of_minutes_until_month[monthCount + 1]:
                month = monthCount + 1
                break
        try:
            day = int((moy - num_of_minutes_until_month[month - 1]) / (60 * 24)) + 1
        except UnboundLocalError:
            raise ValueError(
                "moy must be positive and smaller than 525600. Invalid input %d" % (moy)
            )
        else:
            hour = int((moy / 60) % 24)
            minute = int(moy % 60)

            return cls(month, day, hour, minute, leap_year)

    @classmethod
    def from_date_time_string(cls, datetime_string, leap_year=False):
        """Create Ladybug DateTime from a DateTime string.
        Usage:
            dt = DateTime.from_date_time_string("31 Dec 12:00")
        """
        dt = datetime.strptime(datetime_string, '%d %b %H:%M')
        return cls(dt.month, dt.day, dt.hour, dt.minute, leap_year)

    @property
    def isDateTime(self):
        """Check if data is ladybug data."""
        return True

    @property
    def doy(self):
        """Calculate day of the year for this date time."""
        return self.timetuple().tm_yday

    @property
    def hoy(self):
        """Calculate hour of the year for this date time."""
        return (self.doy - 1) * 24 + self.float_hour

    @property
    def moy(self):
        """Calculate minute of the year for this date time."""
        return self.int_hoy * 60 + self.minute  # minute of the year

    @property
    def float_hour(self):
        """Get hour and minute as a float value (e.g. 6.25 for 6:15)."""
        return self.hour + self.minute / 60.0

    @property
    def int_hoy(self):
        """Calculate hour of the year for this date time as an integer.
        This output assumes the minute is 0.
        """
        return (self.doy - 1) * 24 + self.hour

    @staticmethod
    def _calculate_hour_and_minute(float_hour):
        """Calculate hour and minutes as integers from a float hour."""
        hour, minute = int(float_hour), int(round((float_hour - int(float_hour)) * 60))
        if minute == 60:
            return hour + 1, 0
        else:
            return hour, minute

    def add_minute(self, minute):
        """Create a new DateTime after the minutes are added.
        Args:
            minute: An integer value for minutes.
        """
        _moy = self.moy + int(minute)
        return self.__class__.from_moy(_moy)

    def sub_minute(self, minute):
        """Create a new DateTime after the minutes are subtracted.
        Args:
            minute: An integer value for the number of minutes.
        """
        return self.add_minute(-minute)

    def add_hour(self, hour):
        """Create a new DateTime from this time + timedelta.
        Args:
            hours: A float value in hours (e.g. .5 = half an hour)
        """
        return self.add_minute(hour * 60)

    def sub_hour(self, hour):
        """Create a new DateTime from this time - timedelta.
        Args:
            hour: A float value in hours (e.g. .5 is half an hour and 2 is two hours).
        """
        return self.add_hour(-hour)

    def to_simple_string(self, separator="_"):
        """Return a simplified string."""
        return self.strftime('%d_%b_%H_%M').replace("_", separator)

    def __str__(self):
        """Return date time as a string."""
        return self.strftime('%d %b %H:%M')

    def to_json(self):
        """Get date time as a dictionary."""
        return {'year': self.year,
                'month': self.month,
                'day': self.day,
                'hour': self.hour,
                'minute': self.minute}

    def ToString(self):
        """Overwrite .NET ToString."""
        return self.__str__()

    def __repr__(self):
        """Return date time as a string."""
        return self.__str__()

class Location(object):
    """Ladybug Location.
    Attributes:
        city: Name of the city as a string.
        country: Name of the country as a string.
        latitude: Location latitude between -90 and 90 (Default: 0).
        longitude: Location longitude between -180 (west) and 180 (east) (Default: 0).
        time_zone: Time zone between -12 hours (west) and 12 hours (east) (Default: 0).
        elevation: A number for elevation of the location.
        station_id: Id of the location if the location is represnting a weather station.
        source: Source of data (e.g. TMY, TMY3).
    """

    __slots__ = ("city", "country", "_lat", "_lon", "_tz", "_elev",
                 "station_id", "source")

    def __init__(self, city=None, country=None, latitude=0, longitude=0,
                 time_zone=0, elevation=0, station_id=None, source=None):
        """Create a Ladybug location."""
        self.city = '-' if not city else str(city)
        self.country = '-' if not country else str(country)
        self.latitude = latitude or 0
        self.longitude = longitude or 0
        self.time_zone = time_zone or 0
        self.elevation = elevation or 0
        self.station_id = None if not station_id else str(station_id)
        self.source = source

    @classmethod
    def from_json(cls, loc_json):
        """Create a location from json.
        {
          "city": "-",
          "latitude": 0,
          "longitude": 0,
          "time_zone": 0,
          "elevation": 0
        }
        """
        d = loc_json
        if 'city' not in d:
            d['city'] = None
        if 'country' not in d:
            d['country'] = None
        if 'latitude' not in d:
            d['latitude'] = None
        if 'longitude' not in d:
            d['longitude'] = None
        if 'time_zone' not in d:
            d['time_zone'] = None
        if 'elevation' not in d:
            d['elevation'] = None

        return cls(d['city'], d['country'], d['latitude'], d['longitude'],
                   d['time_zone'], d['elevation'])

    @classmethod
    def from_location(cls, location):
        """Try to create a Ladybug location from a location string.
        Args:
            locationString: Location string
        Usage:
            l = Location.from_location(locationString)
        """
        if not location:
            return cls()
        try:
            if hasattr(location, 'isLocation'):
                # Ladybug location
                return location

            elif hasattr(location, 'Latitude'):
                # Revit's location
                return cls(city=str(location.Name.replace(",", " ")),
                           latitude=location.Latitude,
                           longitude=location.Longitude)

            elif location.startswith('Site:'):
                loc, city, latitude, longitude, time_zone, elevation = \
                    [x.strip() for x in re.findall(r'\r*\n*([^\r\n]*)[,|;]',
                                                   location, re.DOTALL)]
            else:
                try:
                    city, latitude, longitude, time_zone, elevation = \
                        [key.split(":")[-1].strip()
                         for key in location.split(",")]
                except ValueError:
                    # it's just the city name
                    return cls(city=location)

            return cls(city=city, country=None, latitude=latitude,
                       longitude=longitude, time_zone=time_zone,
                       elevation=elevation)

        except Exception as e:
            raise ValueError(
                "Failed to create a Location from %s!\n%s" % (location, e))

    @property
    def isLocation(self):
        """Return Ture."""
        return True

    @property
    def latitude(self):
        """Location latitude."""
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        self._lat = 0 if not lat else float(lat)
        assert -90 <= self._lat <= 90, "latitude should be between -90..90."

    @property
    def longitude(self):
        """Location longitude."""
        return self._lon

    @longitude.setter
    def longitude(self, lon):
        self._lon = 0 if not lon else float(lon)
        assert -180 <= self._lon <= 180, "longitude should be between -180..180."

    @property
    def time_zone(self):
        """Location time zone."""
        return self._tz

    @time_zone.setter
    def time_zone(self, tz):
        self._tz = 0 if not tz else float(tz)
        assert -12 <= self._tz <= 12, "Time zone should be between -12.0..12.0"

    @property
    def elevation(self):
        """Location elevation."""
        return self._elev

    @elevation.setter
    def elevation(self, elev):
        try:
            self._elev = float(elev)
        except TypeError:
            raise ValueError("Failed to convert {} to an elevation".format(elev))

    @property
    def meridian(self):
        """Location meridian west of Greenwich."""
        return -15 * self.time_zone

    def duplicate(self):
        """Duplicate location."""
        return self(self.city, self.country, self.latitude, self.longitude,
                    self.time_zone, self.elevation, self.station_id, self.source)

    @property
    def ep_style_location_string(self):
        """Return EnergyPlus's location string."""
        return "Site:Location,\n  " + \
               self.city + ',\n  ' + \
               str(self.latitude) + ',      !Latitude\n  ' + \
               str(self.longitude) + ',     !Longitude\n  ' + \
               str(self.time_zone) + ',     !Time Zone\n  ' + \
               str(self.elevation) + ';       !Elevation'

    def __str__(self):
        """Return location as a string."""
        return "%s" % (self.ep_style_location_string)

    def ToString(self):
        """Overwrite .NET ToString."""
        return self.__repr__()

    def to_json(self):
        """Create a location from json.
            {
              "city": "-",
              "latitude": 0,
              "longitude": 0,
              "time_zone": 0,
              "elevation": 0
            }
        """
        return {
            "city": self.city,
            "country": self.country,
            "latitude": self.latitude,
            "longitude": self.longitude,
            "time_zone": self.time_zone,
            "elevation": self.elevation
        }

    def __repr__(self):
        """Return location as a string."""
        # Tehran, lat:36, lon:34, tz:3.5, elev:54
        return "%s, lat:%.2f, lon:%.2f, tz:%.1f, elev:%.2f" % (
            self.city, self.latitude, self.longitude,
            self.time_zone, self.elevation)

class Sunpath(object):
    """
    Calculates sun path.
    Attributes:
        latitude: The latitude of the location in degrees. Values must be
            between -90 and 90. Default is set to the equator.
        longitude: The longitude of the location in degrees (Default: 0)
        time_zone: A number representing the time zone of the location you are
            constructing. This can improve the accuracy of the resulting sun
            plot.  The time zone should follow the epw convention and should be
            between -12 and +12, where 0 is at Greenwich, UK, positive values
            are to the East of Greenwich and negative values are to the West.
        north_angle: Angle to north (0-360). 90 is west and 270 is east
            (Default: 0)
        daylight_saving_period: An analysis period for daylight saving.
            (Default: None)
    Usage:
        import ladybug.sunpath as sunpath
        # initiate sunpath
        sp = sunpath.Sunpath(50)
        sun = sp.calculate_sun(1, 1, 12) # calculate sun data for Jan 1 at noon
        print(sun.azimuth, sun.altitude)
    """

    __slots__ = ('_longitude', '_latitude', 'north_angle', 'time_zone',
                 'daylight_saving_period', '_is_leap_year')
    PI = math.pi

    def __init__(self, latitude=0, longitude=0, time_zone=0, north_angle=0,
                 daylight_saving_period=None):
        """Init sunpath.
        Args:
            latitude: The latitude of the location in degrees. Values must be
                between -90 and 90. Default is set to the equator.
            longitude: The longitude of the location in degrees (Default: 0)
            time_zone: A number representing the time zone of the location you are
                constructing. This can improve the accuracy of the resulting sun
                plot.  The time zone should follow the epw convention and should be
                between -12 and +12, where 0 is at Greenwich, UK, positive values
                are to the East of Greenwich and negative values are to the West.
            north_angle: Angle to north (0-360). 90 is west and 270 is east
                (Default: 0).
            daylight_saving_period: An analysis period for daylight saving.
                (Default: None).
        """
        self.time_zone = time_zone
        self.latitude = latitude
        self.longitude = longitude
        self.north_angle = north_angle
        self.daylight_saving_period = daylight_saving_period
        self._is_leap_year = False

    @classmethod
    def from_location(cls, location, north_angle=0, daylight_saving_period=None):
        """Create a sun path from a LBlocation."""
        location = Location.from_location(location)
        return cls(location.latitude, location.longitude,
                   location.time_zone, north_angle, daylight_saving_period)

    @property
    def latitude(self):
        """Get/set latitude in degrees."""
        return math.degrees(self._latitude)

    @latitude.setter
    def latitude(self, value):
        """Set latitude value."""
        self._latitude = math.radians(float(value))
        assert -self.PI / 2 <= self._latitude <= self.PI / 2, \
            "latitude value should be between -90..90."

    @property
    def longitude(self):
        """Get longitude in degrees."""
        return math.degrees(self._longitude)

    @longitude.setter
    def longitude(self, value):
        """Set longitude value in degrees."""
        self._longitude = math.radians(float(value))

        # update time_zone
        if abs((value / 15.0) - self.time_zone) > 1:
            # if time_zone doesn't match the longitude update the time_zone
            self.time_zone = value / 15.0

    @property
    def is_leap_year(self):
        """Indicate is sunpath calculated for a leap year."""
        return self._is_leap_year

    @is_leap_year.setter
    def is_leap_year(self, value):
        """set sunpath to be calculated for a leap year."""
        self._is_leap_year = bool(value)

    def is_daylight_saving_hour(self, datetime):
        """Check if a datetime is a daylight saving time."""
        if not self.daylight_saving_period:
            return False
        return self.daylight_saving_period.isTimeIncluded(datetime.hoy)

    def calculate_sun(self, month, day, hour, is_solar_time=False):
        """Get Sun data for an hour of the year.
        Args:
            month: An integer between 1-12
            day: An integer between 1-31
            hour: A positive number between 0..23
            is_solar_time: A boolean to indicate if the input hour is solar time.
                (Default: False)
        Returns:
            A sun object for this particular time
        """
        datetime = DateTime(month, day, *self._calculate_hour_and_minute(hour),
                            leap_year=self.is_leap_year)
        return self.calculate_sun_from_date_time(datetime, is_solar_time)

    def calculate_sun_from_hoy(self, hoy, is_solar_time=False):
        """Get Sun data for an hour of the year.
        Args:
            datetime: Ladybug datetime
            is_solar_time: A boolean to indicate if the input hour is solar time
                (Default: False).
        Returns:
            A sun object for this particular time
        """
        datetime = DateTime.from_hoy(hoy, self.is_leap_year)
        return self.calculate_sun_from_date_time(datetime, is_solar_time)

    def calculate_sun_from_date_time(self, datetime, is_solar_time=False):
        """Get Sun for an hour of the year.
        This code is originally written by Trygve Wastvedt \
         (Trygve.Wastvedt@gmail.com)
        based on (NOAA) and modified by Chris Mackey and Mostapha Roudsari
        Args:
            datetime: Ladybug datetime
            is_solar_time: A boolean to indicate if the input hour is solar time.
                (Default: False)
        Returns:
            A sun object for this particular time
        """
        # TODO(mostapha): This should be more generic and based on a method
        if datetime.year != 2016 and self.is_leap_year:
            datetime = DateTime(datetime.month, datetime.day, datetime.hour,
                                datetime.minute, True)

        sol_dec, eq_of_time = self._calculate_solar_geometry(datetime)

        hour = datetime.float_hour

        is_daylight_saving = self.is_daylight_saving_hour(datetime.hoy)

        hour = hour + 1 if self.is_daylight_saving_hour(datetime.hoy) else hour

        # minutes
        sol_time = self._calculate_solar_time(hour, eq_of_time, is_solar_time) * 60

        # degrees
        if sol_time / 4 < 0:
            hour_angle = sol_time / 4 + 180
        else:
            hour_angle = sol_time / 4 - 180

        # Degrees
        zenith = math.degrees(math.acos
                              (math.sin(self._latitude) *
                               math.sin(math.radians(sol_dec)) +
                               math.cos(self._latitude) *
                               math.cos(math.radians(sol_dec)) *
                               math.cos(math.radians(hour_angle))))

        altitude = 90 - zenith

        # Approx Atmospheric Refraction
        if altitude > 85:
            atmos_refraction = 0
        else:
            if altitude > 5:
                atmos_refraction = 58.1 / math.tan(math.radians(altitude))

                - 0.07 / (math.tan(math.radians(altitude)))**3
                + 0.000086 / (math.tan(math.radians(altitude)))**5
            else:
                if altitude > -0.575:
                    atmos_refraction = 1735

                    + altitude * (-518.2 + altitude *
                                  (103.4 + altitude *
                                   (-12.79 + altitude * 0.711)))
                else:

                    atmos_refraction = -20.772 / math.tan(
                        math.radians(altitude))

        atmos_refraction /= 3600

        altitude += atmos_refraction

        # Degrees
        if hour_angle > 0:
            azimuth = (math.degrees(
                math.acos(
                    (
                            (math.sin(self._latitude) *
                             math.cos(math.radians(zenith))) -
                            math.sin(math.radians(sol_dec))) /
                    (math.cos(self._latitude) *
                     math.sin(math.radians(zenith)))
                )
            ) + 180) % 360
        else:
            azimuth = (540 - math.degrees(math.acos((
                                                            (math.sin(self._latitude) *
                                                             math.cos(math.radians(zenith))) -
                                                            math.sin(math.radians(sol_dec))) /
                                                    (math.cos(self._latitude) *
                                                     math.sin(math.radians(zenith))))
                                          )) % 360

        altitude = math.radians(altitude)
        azimuth = math.radians(azimuth)
        # create the sun for this hour
        return Sun(datetime, altitude, azimuth, is_solar_time, is_daylight_saving,
                   self.north_angle)

    def calculate_sunrise_sunset(self, month, day, depression=0.833,
                                 is_solar_time=False):
        """Calculate sunrise, noon and sunset.
        Return:
            A dictionary. Keys are ("sunrise", "noon", "sunset")
        """
        datetime = DateTime(month, day, hour=12, leap_year=self.is_leap_year)

        return self.calculate_sunrise_sunset_from_datetime(datetime,
                                                           depression,
                                                           is_solar_time)

    # TODO: implement solar time
    def calculate_sunrise_sunset_from_datetime(self, datetime, depression=0.833,
                                               is_solar_time=False):
        """Calculate sunrise, sunset and noon for a day of year."""
        # TODO(mostapha): This should be more generic and based on a method
        if datetime.year != 2016 and self.is_leap_year:
            datetime = DateTime(datetime.month, datetime.day, datetime.hour,
                                datetime.minute, True)
        sol_dec, eq_of_time = self._calculate_solar_geometry(datetime)
        # calculate sunrise and sunset hour
        if is_solar_time:
            noon = .5
        else:
            noon = (720 -
                    4 * math.degrees(self._longitude) -
                    eq_of_time +
                    self.time_zone * 60
                    ) / 1440.0

        try:
            sunrise_hour_angle = self._calculate_sunrise_hour_angle(
                sol_dec, depression)
        except ValueError:
            # no sun rise and sunset at this hour
            noon = 24 * noon
            return {
                "sunrise": None,
                "noon": DateTime(datetime.month, datetime.day,
                                 *self._calculate_hour_and_minute(noon),
                                 leap_year=self.is_leap_year),
                "sunset": None
            }
        else:
            sunrise = noon - sunrise_hour_angle * 4 / 1440.0
            sunset = noon + sunrise_hour_angle * 4 / 1440.0
            noon = 24 * noon
            sunrise = 24 * sunrise
            sunset = 24 * sunset

            return {
                "sunrise": DateTime(datetime.month, datetime.day,
                                    *self._calculate_hour_and_minute(sunrise),
                                    leap_year=self.is_leap_year),
                "noon": DateTime(datetime.month, datetime.day,
                                 *self._calculate_hour_and_minute(noon),
                                 leap_year=self.is_leap_year),
                "sunset": DateTime(datetime.month, datetime.day,
                                   *self._calculate_hour_and_minute(sunset),
                                   leap_year=self.is_leap_year)
            }

    def _calculate_solar_geometry(self, datetime):
        """Calculate Solar geometry for an hour of the year.
        Attributes:
            datetime: A Ladybug datetime
        Returns:
            Solar declination: Solar declination in radians
            eq_of_time: Equation of time as minutes
        """
        month = datetime.month
        day = datetime.day
        hour = datetime.hour
        minute = datetime.minute
        year = 2016 if self.is_leap_year else 2017

        def find_fraction_of_24(hour, minute):
            """
            This function calculates the fraction of the 24 hour
            the provided time represents
            1440 is total the number of minutes in a 24 hour cycle.
            args
                hour: Integer. Hour between 0 - 23
                minute: Integer. Minute between 0 - 59
            return: Float.
                The fraction of the 24 hours the provided time represents
            """
            return round((minute + hour * 60) / 1440.0, 2)

        def days_from_010119(year, month, day):
            """
            This function calculates the number of days from 01-01-1900 \
            to the provided date
            args :
                year: Integer. The year in the date
                month: Integer. The month in the date
                day: Integer. The date
            return: The number of days from 01-01-1900 to the date provided
            """

            # Making a list of years from the year 1900
            years = range(1900, year)

            def is_leap_year(year):
                """Determine whether a year is a leap year."""
                return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

            # Number of days in a year are 366 if it is a leap year
            days_in_year = []
            for item in years:
                if is_leap_year(item):
                    days_in_year.append(366)
                else:
                    days_in_year.append(365)

            # Making the total of all the days in preceding years
            days_in_precending_years = 0
            for days in days_in_year:
                days_in_precending_years += days

            if is_leap_year(year):
                month_dict = {1: 31, 2: 29, 3: 31, 4: 30, 5: 31, 6: 30,
                              7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
            else:
                month_dict = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30,
                              7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}

            """Making the total of all the days in preceding months\
            in the same year"""
            keys = tuple(month_dict.keys())
            days_in_precending_months = 0
            for i in range(month - 1):
                days_in_precending_months += month_dict[keys[i]]

            return days_in_precending_years + days_in_precending_months + day + 1

        julian_day = days_from_010119(year, month, day) + 2415018.5 + \
                     find_fraction_of_24(hour, minute) - (float(self.time_zone) / 24)

        julian_century = (julian_day - 2451545) / 36525

        # degrees
        geom_mean_long_sun = (280.46646 + julian_century *
                              (36000.76983 + julian_century * 0.0003032)
                              ) % 360
        # degrees
        geom_mean_anom_sun = 357.52911 + julian_century * \
                             (35999.05029 - 0.0001537 * julian_century)

        eccent_orbit = 0.016708634 - julian_century * \
                       (0.000042037 + 0.0000001267 * julian_century)

        sun_eq_of_ctr = math.sin(
            math.radians(geom_mean_anom_sun)) * \
                        (1.914602 - julian_century * (0.004817 + 0.000014 * julian_century)
                         ) + \
                        math.sin(math.radians(2 * geom_mean_anom_sun)) * \
                        (0.019993 - 0.000101 * julian_century) + \
                        math.sin(math.radians(3 * geom_mean_anom_sun)) * \
                        0.000289

        # degrees
        sun_true_long = geom_mean_long_sun + sun_eq_of_ctr

        # degrees
        sun_app_long = sun_true_long - 0.00569 - 0.00478 * \
                       math.sin(math.radians(125.04 - 1934.136 * julian_century))

        # degrees
        mean_obliq_ecliptic = 23 + \
                              (26 + ((21.448 - julian_century * (46.815 + julian_century *
                                                                 (0.00059 - julian_century *
                                                                  0.001813)))) / 60) / 60

        # degrees
        oblique_corr = mean_obliq_ecliptic + 0.00256 * \
                       math.cos(math.radians(125.04 - 1934.136 * julian_century))

        # RADIANS
        sol_dec = math.degrees(math.asin(math.sin(math.radians(oblique_corr)) *
                                         math.sin(math.radians(sun_app_long))))

        var_y = math.tan(math.radians(oblique_corr / 2)) * \
                math.tan(math.radians(oblique_corr / 2))

        # minutes
        eq_of_time = 4 \
                     * math.degrees(
            var_y * math.sin(2 * math.radians(geom_mean_long_sun)) -
            2 * eccent_orbit * math.sin(math.radians(geom_mean_anom_sun)) +
            4 * eccent_orbit * var_y *
            math.sin(math.radians(geom_mean_anom_sun)) *
            math.cos(2 * math.radians(geom_mean_long_sun)) -
            0.5 * (var_y ** 2) *
            math.sin(4 * math.radians(geom_mean_long_sun)) -
            1.25 * (eccent_orbit ** 2) *
            math.sin(2 * math.radians(geom_mean_anom_sun))
        )

        return sol_dec, eq_of_time

    def _calculate_sunrise_hour_angle(self, solar_dec, depression=0.833):
        """Calculate hour angle for sunrise time in degrees."""

        hour_angle_arg = math.degrees(math.acos(
            math.cos(math.radians(90 + depression)) /
            (math.cos(math.radians(self.latitude)) * math.cos(
                math.radians(solar_dec))) -
            math.tan(math.radians(self.latitude)) *
            math.tan(math.radians(solar_dec))
        ))

        return hour_angle_arg

    def _calculate_solar_time(self, hour, eq_of_time, is_solar_time):
        """Calculate Solar time for an hour."""
        if is_solar_time:
            return hour

        return (
                       (hour * 60 + eq_of_time + 4 * math.degrees(self._longitude) -
                        60 * self.time_zone) % 1440) / 60

    def _calculate_solar_time_by_doy(self, hour, doy):
        """This is how radiance calculates solar time.
        This is a place holder and \
        need to be validated against calculateSolarTime.
        """
        raise NotImplementedError()
        return (0.170 * math.sin((4 * math.pi / 373) * (doy - 80)) -
                0.129 * math.sin((2 * math.pi / 355) * (doy - 8)) +
                12 * (-(15 * self.time_zone) - self.longitude) / math.pi)

    @staticmethod
    def _calculate_hour_and_minute(float_hour):
        """Calculate hour and minutes as integers from a float hour."""
        hour = int(float_hour)
        minute = int(round((float_hour - int(float_hour)) * 60))

        if minute == 60:
            return hour + 1, 0
        else:
            return hour, minute

    def draw_sunpath(self,
                     hoys=None,
                     origin=None,
                     scale=1, sun_scale=1, annual=True, rem_night=True):
        """Create sunpath geometry. \
        This method should only be used from the + libraries.
        Args:
            hoys: An optional list of hours of the year(default: None).
            origin: Sunpath origin(default: (0, 0, 0)).
            scale: Sunpath scale(default: 1).
            sun_scale: Scale for the sun spheres(default: 1).
            annual: Set to True to draw an annual sunpath.
                Otherwise a daily sunpath is drawn.
            rem_night: Remove suns which are under the horizon(night!).
        Returns:
            base_curves: A collection of curves for base plot.
            analemma_curves: A collection of analemma_curves.
            daily_curves: A collection of daily_curves.
            suns: A list of suns.
        """
        # check and make sure the call is coming from inside a plus library
        assert ladybug.isplus, \
            '"draw_sunpath" method can only be used in the [+] libraries.'
        hoys = hoys or ()
        origin = origin or (0, 0, 0)
        try:
            origin = tuple(origin)
        except TypeError as e:
            # dynamo
            try:
                origin = origin.X, origin.Y, origin.Z
            except AttributeError:
                raise TypeError(str(e))

        scale = scale or 1
        sun_scale = sun_scale or 1
        assert annual or hoys, 'For daily sunpath you need to provide at least one hour.'

        radius = 200 * scale

        # draw base circles and lines
        base_curves = plus.base_curves(origin, radius, self.north_angle)
        # draw analemma
        # calculate date times for analemma curves
        if annual:
            asuns = self._analemma_suns()
            analemma_curves = plus.analemma_curves(asuns, origin, radius)
        else:
            analemma_curves = ()

        # add sun spheres
        if hoys:
            suns = tuple(self.calculate_sun_from_hoy(hour) for hour in hoys)
        else:
            suns = ()

        if rem_night:
            suns = tuple(sun for sun in suns if sun.is_during_day)

        sun_geos = plus.sun_geometry(suns, origin, radius)

        # draw daily sunpath
        if annual:
            dts = (DateTime(m, 21) for m in xrange(1, 13))
        else:
            dts = (sun.datetime for sun in suns)

        dsuns = self._daily_suns(dts)
        daily_curves = plus.daily_curves(dsuns, origin, radius)

        SPGeo = namedtuple(
            'SunpathGeo',
            ('compass_curves',
             'analemma_curves',
             'daily_curves',
             'suns',
             'sun_geos'))

        # return outputs
        return SPGeo(base_curves, analemma_curves, daily_curves, suns, sun_geos)

    def _analemma_position(self, hour):
        """Check what the analemma position is for an hour.
        This is useful for calculating hours of analemma curves.
        Returns:
            -1 if always night,
            0 if both day and night,
            1 if always day.
        """
        # check for 21 dec and 21 jun
        low = self.calculate_sun(12, 21, hour).is_during_day
        high = self.calculate_sun(6, 21, hour).is_during_day

        if low and high:
            return 1
        elif low or high:
            return 0
        else:
            return -1

    def _analemma_suns(self):
        """Calculate times that should be used for drawing analemma_curves.
        Returns:
            A list of list of analemma suns.
        """
        for h in xrange(0, 24):
            if self._analemma_position(h) < 0:
                continue
            elif self._analemma_position(h) == 0:
                chours = []
                # this is an hour that not all the hours are day or night
                prevhour = self.latitude <= 0
                num_of_days = 8760 if not self.is_leap_year else 8760 + 24
                for hoy in xrange(h, num_of_days, 24):
                    thishour = self.calculate_sun_from_hoy(hoy).is_during_day
                    if thishour != prevhour:
                        if not thishour:
                            hoy -= 24
                        dt = DateTime.from_hoy(hoy, self.is_leap_year)
                        chours.append((dt.month, dt.day, dt.hour))
                    prevhour = thishour
                tt = []
                for hcount in range(int(len(chours) / 2)):
                    st = chours[2 * hcount]
                    en = chours[2 * hcount + 1]
                    if self.latitude >= 0:
                        tt = [self.calculate_sun(*st)] + \
                             [self.calculate_sun(st[0], d, h)
                              for d in xrange(st[1] + 1, 29, 7)] + \
                             [self.calculate_sun(m, d, h)
                              for m in xrange(st[0] + 1, en[0])
                              for d in xrange(3, 29, 7)] + \
                             [self.calculate_sun(en[0], d, h)
                              for d in xrange(3, en[1], 7)] + \
                             [self.calculate_sun(*en)]
                    else:
                        tt = [self.calculate_sun(*en)] + \
                             [self.calculate_sun(en[0], d, h)
                              for d in xrange(en[1] + 1, 29, 7)] + \
                             [self.calculate_sun(m, d, h) for m in xrange(en[0] +
                                                                          1, 13)
                              for d in xrange(3, 29, 7)] + \
                             [self.calculate_sun(m, d, h) for m in xrange(1, st[0])
                              for d in xrange(3, 29, 7)] + \
                             [self.calculate_sun(st[0], d, h)
                              for d in xrange(3, st[1], 7)] + \
                             [self.calculate_sun(*st)]
                    yield tt
            else:
                yield tuple(self.calculate_sun((m % 12) + 1, d, h)
                            for m in xrange(0, 13) for d in (7, 14, 21))[:-2]

    def _daily_suns(self, datetimes):
        """Get sun curve for multiple days of the year."""
        for dt in datetimes:
            # calculate sunrise sunset and noon
            nss = self.calculate_sunrise_sunset(dt.month, dt.day)
            dts = tuple(nss[k] for k in ('sunrise', 'noon', 'sunset'))
            if dts[0] is None:
                # circle
                yield (self.calculate_sun(dt.month, dt.day, h) for h in (0, 12,
                                                                         15)), \
                      False
            else:
                # Arc
                yield (self.calculate_sun_from_date_time(dt) for dt in dts), True

class Sun(object):
    """Sun.
    Attributes:
        datetime: A DateTime that represents the datetime for this sun_vector
        altitude: Solar Altitude in **radians**
        azimuth: Solar Azimuth in **radians**
        is_solar_time: A Boolean that indicates if datetime represents the solar
            time.
        is_daylight_saving: A Boolean that indicates if datetime is calculated
            for Daylight saving period
        north_angle: North angle of the sunpath in Degrees. This will be only
            used to calculate the solar vector.
    """

    __slots__ = ('_datetime', '_altitude', '_azimuth', '_is_solar_time',
                 '_is_daylight_saving', '_north_angle', '_hourlyData', '_data',
                 '_sun_vector')
    PI = math.pi

    def __init__(self, datetime, altitude, azimuth, is_solar_time,
                 is_daylight_saving, north_angle, data=None):
        """Init sun."""
        assert hasattr(datetime, 'isDateTime'), \
            "datetime must be a DateTime (not {})".format(type(datetime))
        self._datetime = datetime  # read-only

        assert -self.PI <= altitude <= self.PI, \
            "altitude({}) must be between {} and {}." \
                .format(altitude, -self.PI, self.PI)

        self._altitude = altitude  # read-only

        assert -2 * self.PI <= azimuth <= 2 * self.PI, \
            "azimuth({}) should be between {} and {}." \
                .format(azimuth, -self.PI, self.PI)

        self._azimuth = azimuth  # read-only
        self._is_solar_time = is_solar_time
        self._is_daylight_saving = is_daylight_saving
        # useful to calculate sun vector - sun angle is in degrees
        self._north_angle = north_angle
        self.data = data  # Place holder for hourly data

        self._calculate_sun_vector()

    @property
    def datetime(self):
        """Return datetime."""
        return self._datetime

    @property
    def north_angle(self):
        """Return north angle for +YAxis."""
        return self._north_angle

    @property
    def hoy(self):
        """Return Hour of the year."""
        return self._datetime.hoy

    @property
    def altitude(self):
        """Return solar altitude in degrees."""
        return math.degrees(self._altitude)

    @property
    def azimuth(self):
        """Return solar azimuth in degrees."""
        return math.degrees(self._azimuth)

    @property
    def altitude_in_radians(self):
        """Return solar altitude in radians."""
        return self._altitude

    @property
    def azimuth_in_radians(self):
        """Return solar azimuth in radians."""
        return self._azimuth

    @property
    def is_solar_time(self):
        """Return a Boolean that indicates is datetime is solar time."""
        return self._is_solar_time

    @property
    def is_daylight_saving(self):
        """Return a Boolean that indicates is datetime is solar time."""
        return self._is_daylight_saving

    @property
    def data(self):
        """Get or set data to this sun position."""
        return self._data

    @data.setter
    def data(self, d):
        self._data = d

    @property
    def is_during_day(self):
        """Check if this sun position is during day."""
        # sun vector is flipped to look to the center
        return self.sun_vector.z <= 0

    @property
    def sun_vector(self):
        """Sun vector for this sun.
        Sun vector faces downward(e.g. z will be negative.)
        """
        return self._sun_vector

    def _calculate_sun_vector(self):
        """Calculate sun vector for this sun."""
        z_axis = Vector3(0., 0., -1.)
        x_axis = Vector3(1., 0., 0.)
        north_vector = Vector3(0., 1., 0.)

        # rotate north vector based on azimuth, altitude, and north
        _sun_vector = north_vector \
            .rotate_around(x_axis, self.altitude_in_radians) \
            .rotate_around(z_axis, self.azimuth_in_radians) \
            .rotate_around(z_axis, math.radians(-1 * self.north_angle))

        _sun_vector.normalize()
        try:
            _sun_vector.flip()
        except AttributeError:
            # euclid3
            _sun_vector = Vector3(-1 * _sun_vector.x,
                                  -1 * _sun_vector.y,
                                  -1 * _sun_vector.z)

        self._sun_vector = _sun_vector

    def ToString(self):
        """Overwrite .NET ToString method."""
        return self.__repr__()

    def __repr__(self):
        """Sun representation."""
        # sun at datetime (X, Y, Z)
        return "Sun at {} (x:{}, y:{}, z:{})".format(
            self.datetime,
            self.sun_vector.x,
            self.sun_vector.y,
            self.sun_vector.z
        )


# Within this comment is some code partially stolen from PySolar
####################################################################################################

# class Sun:
#     """Get the location of the sun from a date, time and location"""
#
#     from datetime import datetime, timedelta, timezone
#     import math
#
#
#     def __init__(self, latitude: float = 0, longitude: float = 0, elevation: float = 0, year: int = 1989,
#                  month: int = 12, day: int = 13, hour: int = 14, minute: int = 3,
#                  second: int = 0, microsecond: int = 0, utc_offset: float = 0):
#         """
#         :param latitude:
#         :param longitude:
#         :param elevation:
#         :param year:
#         :param month:
#         :param day:
#         :param hour:
#         :param second:
#         :param microsecond:
#         :param utc_offset:
#         """
#
#         self.latitude = latitude
#         self.longitude = longitude
#         self.elevation = elevation
#         self.year = year
#         self.month = month
#         self.day = day
#         self.hour = hour
#         self.minute = minute
#         self.second = second
#         self.microsecond = microsecond
#         self.utc_offset = utc_offset
#
#         self.dt = self.datetime(year, month, day, hour, minute, second, microsecond, tzinfo=self.timezone(self.timedelta(hours=utc_offset)))
#         self.dtt = self.dt.utctimetuple()
#         self.yday = self.dtt.tm_yday
#
#         self.earth_axis_inclination = 23.45  # degrees
#
#         # Constants
#         self.aberration_coefficients = {'ArgumentOfLatitudeOfMoon': [93.27191, 483202.017538, -0.0036825, 327270.0],
#                                         'LongitudeOfAscendingNode': [125.04452, -1934.136261, 0.0020708, 450000.0],
#                                         'MeanElongationOfMoon': [297.85036, 445267.111480, -0.0019142, 189474.0],
#                                         'MeanAnomalyOfMoon': [134.96298, 477198.867398, 0.0086972, 56250.0],
#                                         'MeanAnomalyOfSun': [357.52772, 35999.050340, -0.0001603, -300000.0]}
#         self.earth_radius = 6378140.0  # meters
#         self.seconds_per_day = 86400
#         self.standard_pressure = 101325.00  # pascals
#         self.standard_temperature = 288.15  # kelvin
#         self.celsius_offset = 273.15  # subtract from kelvin to get deg C, add to deg C to get kelvin
#         self.earth_temperature_lapse_rate = -0.0065  # change in temperature with height, kelvin/metre
#         self.air_gas_constant = 8.31432  # N*m/s^2
#         self.earth_gravity = 9.80665  # m/s^2 or N/kg
#         self.earth_atmosphere_molar_mass = 0.0289644  # kg/mol
#         self.aberration_sin_terms = [(0, 0, 0, 0, 1), (-2, 0, 0, 2, 2), (0, 0, 0, 2, 2), (0, 0, 0, 0, 2), (0, 1, 0, 0, 0),
#                                      (0, 0, 1, 0, 0), (-2, 1, 0, 2, 2), (0, 0, 0, 2, 1), (0, 0, 1, 2, 2), (-2, -1, 0, 2, 2),
#                                      (-2, 0, 1, 0, 0), (-2, 0, 0, 2, 1), (0, 0, -1, 2, 2), (2, 0, 0, 0, 0), (0, 0, 1, 0, 1),
#                                      (2, 0, -1, 2, 2), (0, 0, -1, 0, 1), (0, 0, 1, 2, 1), (-2, 0, 2, 0, 0), (0, 0, -2, 2, 1),
#                                      (2, 0, 0, 2, 2), (0, 0, 2, 2, 2), (0, 0, 2, 0, 0), (-2, 0, 1, 2, 2), (0, 0, 0, 2, 0),
#                                      (-2, 0, 0, 2, 0), (0, 0, -1, 2, 1), (0, 2, 0, 0, 0), (2, 0, -1, 0, 1), (-2, 2, 0, 2, 2),
#                                      (0, 1, 0, 0, 1), (-2, 0, 1, 0, 1), (0, -1, 0, 0, 1), (0, 0, 2, -2, 0), (2, 0, -1, 2, 1),
#                                      (2, 0, 1, 2, 2), (0, 1, 0, 2, 2), (-2, 1, 1, 0, 0), (0, -1, 0, 2, 2), (2, 0, 0, 2, 1),
#                                      (2, 0, 1, 0, 0), (-2, 0, 2, 2, 2), (-2, 0, 1, 2, 1), (2, 0, -2, 0, 1), (2, 0, 0, 0, 1),
#                                      (0, -1, 1, 0, 0), (-2, -1, 0, 2, 1), (-2, 0, 0, 0, 1), (0, 0, 2, 2, 1), (-2, 0, 2, 0, 1),
#                                      (-2, 1, 0, 2, 1), (0, 0, 1, -2, 0), (-1, 0, 1, 0, 0), (-2, 1, 0, 0, 0), (1, 0, 0, 0, 0),
#                                      (0, 0, 1, 2, 0), (0, 0, -2, 2, 2), (-1, -1, 1, 0, 0), (0, 1, 1, 0, 0), (0, -1, 1, 2, 2),
#                                      (2, -1, -1, 2, 2), (0, 0, 3, 2, 2), (2, -1, 0, 2, 2)]
#         self.nutation_coefficients = [(-171996, -174.2, 92025, 8.9), (-13187, -1.6, 5736, -3.1), (-2274, -0.2, 977, -0.5),
#                                       (2062, 0.2, -895, 0.5), (1426, -3.4, 54, -0.1), (712, 0.1, -7, 0), (-517, 1.2, 224, -0.6),
#                                       (-386, -0.4, 200, 0), (-301, 0, 129, -0.1), (217, -0.5, -95, 0.3), (-158, 0, 0, 0),
#                                       (129, 0.1, -70, 0), (123, 0, -53, 0), (63, 0, 0, 0), (63, 0.1, -33, 0), (-59, 0, 26, 0),
#                                       (-58, -0.1, 32, 0), (-51, 0, 27, 0), (48, 0, 0, 0), (46, 0, -24, 0), (-38, 0, 16, 0),
#                                       (-31, 0, 13, 0), (29, 0, 0, 0), (29, 0, -12, 0), (26, 0, 0, 0), (-22, 0, 0, 0),
#                                       (21, 0, -10, 0), (17, -0.1, 0, 0), (16, 0, -8, 0), (-16, 0.1, 7, 0), (-15, 0, 9, 0),
#                                       (-13, 0, 7, 0), (-12, 0, 6, 0), (11, 0, 0, 0), (-10, 0, 5, 0), (-8, 0, 3, 0), (7, 0, -3, 0),
#                                       (-7, 0, 0, 0), (-7, 0, 3, 0), (-7, 0, 3, 0), (6, 0, 0, 0), (6, 0, -3, 0), (6, 0, -3, 0),
#                                       (-6, 0, 3, 0), (-6, 0, 3, 0), (5, 0, 0, 0), (-5, 0, 3, 0), (-5, 0, 3, 0), (-5, 0, 3, 0),
#                                       (4, 0, 0, 0), (4, 0, 0, 0), (4, 0, 0, 0), (-4, 0, 0, 0), (-4, 0, 0, 0), (-4, 0, 0, 0),
#                                       (3, 0, 0, 0), (-3, 0, 0, 0), (-3, 0, 0, 0), (-3, 0, 0, 0), (-3, 0, 0, 0), (-3, 0, 0, 0),
#                                       (-3, 0, 0, 0), (-3, 0, 0, 0)]
#         self.heliocentric_longitude_coefficients = [
#             [(175347046.0, 0, 0), (3341656.0, 4.6692568, 6283.07585), (34894.0, 4.6261, 12566.1517),
#              (3497.0, 2.7441, 5753.3849), (3418.0, 2.8289, 3.5231), (3136.0, 3.6277, 77713.7715), (2676.0, 4.4181, 7860.4194),
#              (2343.0, 6.1352, 3930.2097), (1324.0, 0.7425, 11506.7698), (1273.0, 2.0371, 529.691), (1199.0, 1.1096, 1577.3435),
#              (990, 5.233, 5884.927), (902, 2.045, 26.298), (857, 3.508, 398.149), (780, 1.179, 5223.694),
#              (753, 2.533, 5507.553), (505, 4.583, 18849.228), (492, 4.205, 775.523), (357, 2.92, 0.067),
#              (317, 5.849, 11790.629), (284, 1.899, 796.298), (271, 0.315, 10977.079), (243, 0.345, 5486.778),
#              (206, 4.806, 2544.314), (205, 1.869, 5573.143), (202, 2.4458, 6069.777), (156, 0.833, 213.299),
#              (132, 3.411, 2942.463), (126, 1.083, 20.775), (115, 0.645, 0.98), (103, 0.636, 4694.003), (102, 0.976, 15720.839),
#              (102, 4.267, 7.114), (99, 6.21, 2146.17), (98, 0.68, 155.42), (86, 5.98, 161000.69), (85, 1.3, 6275.96),
#              (85, 3.67, 71430.7), (80, 1.81, 17260.15), (79, 3.04, 12036.46), (71, 1.76, 5088.63), (74, 3.5, 3154.69),
#              (74, 4.68, 801.82), (70, 0.83, 9437.76), (62, 3.98, 8827.39), (61, 1.82, 7084.9), (57, 2.78, 6286.6),
#              (56, 4.39, 14143.5), (56, 3.47, 6279.55), (52, 0.19, 12139.55), (52, 1.33, 1748.02), (51, 0.28, 5856.48),
#              (49, 0.49, 1194.45), (41, 5.37, 8429.24), (41, 2.4, 19651.05), (39, 6.17, 10447.39), (37, 6.04, 10213.29),
#              (37, 2.57, 1059.38), (36, 1.71, 2352.87), (36, 1.78, 6812.77), (33, 0.59, 17789.85), (30, 0.44, 83996.85),
#              (30, 2.74, 1349.87), (25, 3.16, 4690.48)],
#             [(628331966747.0, 0, 0), (206059.0, 2.678235, 6283.07585), (4303.0, 2.6351, 12566.1517), (425.0, 1.59, 3.523),
#              (119.0, 5.796, 26.298), (109.0, 2.966, 1577.344), (93, 2.59, 18849.23), (72, 1.14, 529.69), (68, 1.87, 398.15),
#              (67, 4.41, 5507.55), (59, 2.89, 5223.69), (56, 2.17, 155.42), (45, 0.4, 796.3), (36, 0.47, 775.52),
#              (29, 2.65, 7.11), (21, 5.34, 0.98), (19, 1.85, 5486.78), (19, 4.97, 213.3), (17, 2.99, 6275.96),
#              (16, 0.03, 2544.31), (16, 1.43, 2146.17), (15, 1.21, 10977.08), (12, 2.83, 1748.02), (12, 3.26, 5088.63),
#              (12, 5.27, 1194.45), (12, 2.08, 4694), (11, 0.77, 553.57), (10, 1.3, 3286.6), (10, 4.24, 1349.87),
#              (9, 2.7, 242.73), (9, 5.64, 951.72), (8, 5.3, 2352.87), (6, 2.65, 9437.76), (6, 4.67, 4690.48)],
#             [(52919.0, 0, 0), (8720.0, 1.0721, 6283.0758), (309.0, 0.867, 12566.152), (27, 0.05, 3.52), (16, 5.19, 26.3),
#              (16, 3.68, 155.42), (10, 0.76, 18849.23), (9, 2.06, 77713.77), (7, 0.83, 775.52), (5, 4.66, 1577.34),
#              (4, 1.03, 7.11), (4, 3.44, 5573.14), (3, 5.14, 796.3), (3, 6.05, 5507.55), (3, 1.19, 242.73), (3, 6.12, 529.69),
#              (3, 0.31, 398.15), (3, 2.28, 553.57), (2, 4.38, 5223.69), (2, 3.75, 0.98)],
#             [(289.0, 5.844, 6283.076), (35, 0, 0), (17, 5.49, 12566.15), (3, 5.2, 155.42), (1, 4.72, 3.52), (1, 5.3, 18849.23),
#              (1, 5.97, 242.73)], [(114.0, 3.142, 0), (8, 4.13, 6283.08), (1, 3.84, 12566.15)], [(1, 3.14, 0)], ]
#         self.heliocentric_latitude_coefficients = [
#             [(280.0, 3.199, 84334.662), (102.0, 5.422, 5507.553), (80, 3.88, 5223.69), (44, 3.7, 2352.87), (32, 4, 1577.34)],
#             [(9, 3.9, 5507.55), (6, 1.73, 5223.69)], ]
#         self.sun_earth_distance_coefficients = [
#             [(100013989.0, 0, 0), (1670700.0, 3.0984635, 6283.07585), (13956.0, 3.05525, 12566.1517),
#              (3084.0, 5.1985, 77713.7715), (1628.0, 1.1739, 5753.3849), (1576.0, 2.8469, 7860.4194), (925.0, 5.453, 11506.77),
#              (542.0, 4.564, 3930.21), (472.0, 3.661, 5884.927), (346.0, 0.964, 5507.553), (329.0, 5.9, 5223.694),
#              (307.0, 0.299, 5573.143), (243.0, 4.273, 11790.629), (212.0, 5.847, 1577.344), (186.0, 5.022, 10977.079),
#              (175.0, 3.012, 18849.228), (110.0, 5.055, 5486.778), (98, 0.89, 6069.78), (86, 5.69, 15720.84),
#              (86, 1.27, 161000.69), (85, 0.27, 17260.15), (63, 0.92, 529.69), (57, 2.01, 83996.85), (56, 5.24, 71430.7),
#              (49, 3.25, 2544.31), (47, 2.58, 775.52), (45, 5.54, 9437.76), (43, 6.01, 6275.96), (39, 5.36, 4694),
#              (38, 2.39, 8827.39), (37, 0.83, 19651.05), (37, 4.9, 12139.55), (36, 1.67, 12036.46), (35, 1.84, 2942.46),
#              (33, 0.24, 7084.9), (32, 0.18, 5088.63), (32, 1.78, 398.15), (28, 1.21, 6286.6), (28, 1.9, 6279.55),
#              (26, 4.59, 10447.39)],
#             [(103019.0, 1.10749, 6283.07585), (1721.0, 1.0644, 12566.1517), (702.0, 3.142, 0), (32, 1.02, 18849.23),
#              (31, 2.84, 5507.55), (25, 1.32, 5223.69), (18, 1.42, 1577.34), (10, 5.91, 10977.08), (9, 1.42, 6275.96),
#              (9, 0.27, 5486.78)],
#             [(4359.0, 5.7846, 6283.0758), (124.0, 5.579, 12566.152), (12, 3.14, 0), (9, 3.63, 77713.77), (6, 1.87, 5573.14),
#              (3, 5.47, 18849)], [(145.0, 4.273, 6283.076), (7, 3.92, 12566.15), ], [(4, 2.56, 6283.08)], ]
#         self.julian_day_offset = 1721425 - 0.5  # add to datetime.datetime.toordinal() to get Julian day number
#         self.gregorian_day_offset = 719163  # number of days to add to datetime.datetime.timestamp() / seconds_per_day to agree with datetime.datetime.toordinal()
#         self.tt_offset = 32.184  # seconds to add to TAI to get TT
#         self.leap_seconds_base_year = 1972
#         self.leap_seconds_adjustments = [
#             # two entries per year starting from 1972, first for 23:59:59 June 30,
#             # second for 23:59:59 December 31. +1 indicates that 23:59:60 follows,
#             # -1 indicates that 23:59:59 does not exist, not that the latter has ever occurred.
#             # source: https://www.nist.gov/pml/time-and-frequency-division/atomic-standards/leap-second-and-ut1-utc-information
#             (+1, +1), (0, +1), (0, +1), (0, +1), (0, +1), (0, +1), (0, +1), (0, +1), (0, 0), (+1, 0), (+1, 0), (+1, 0), (0, 0),
#             (+1, 0), (0, 0), (0, +1), (0, 0), (0, +1), (0, +1), (0, 0), (+1, 0), (+1, 0), (+1, 0), (0, +1), (0, 0), (+1, 0),
#             (+1, 0), (0, 0), (0, +1), (0, 0), (0, +1), (0, +1), (0, 0), (+1, 0), (+1, 0), (+1, 0), (0, +1), (0, 0), (+1, 0),
#             (0, +1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, +1), (0, 0), (0, 0), (0, +1), (0, 0), (0, 0), (0, 0),
#             (+1, 0), (0, 0), (0, 0), (+1, 0), (0, +1), (0, 0), (0, 0),
#         ]
#
#         # table of values to add to UT1 to get TT (to date), generated by util/get_delta_t script
#         self.delta_t_base_year = 1973
#         self.delta_t_base_month = 2
#         self.delta_t = [
#             [43.4724, 43.5648, 43.6737, 43.7782, 43.8763, 43.9562, 44.0315, 44.1132, 44.1982, 44.2952, 44.3936],
#             [44.4841, 44.5646, 44.6425, 44.7386, 44.8370, 44.9302, 44.9986, 45.0584, 45.1284, 45.2064, 45.2980,
#              45.3897],
#             [45.4761, 45.5633, 45.6450, 45.7375, 45.8284, 45.9133, 45.9820, 46.0408, 46.1067, 46.1825,
#              46.2789, 46.3713],
#             [46.4567, 46.5445, 46.6311, 46.7302, 46.8284, 46.9247, 46.9970, 47.0709, 47.1451, 47.2362,
#              47.3413, 47.4319],
#             [47.5214, 47.6049, 47.6837, 47.7781, 47.8771, 47.9687, 48.0348, 48.0942, 48.1608, 48.2460,
#              48.3439, 48.4355],
#             [48.5344, 48.6325, 48.7294, 48.8365, 48.9353, 49.0319, 49.1013, 49.1591, 49.2286, 49.3070,
#              49.4018, 49.4945],
#             [49.5862, 49.6805, 49.7602, 49.8556, 49.9489, 50.0347, 50.1019, 50.1622, 50.2260, 50.2968,
#              50.3831, 50.4599],
#             [50.5387, 50.6161, 50.6866, 50.7658, 50.8454, 50.9187, 50.9761, 51.0278, 51.0843, 51.1538,
#              51.2319, 51.3063],
#             [51.3808, 51.4526, 51.5160, 51.5985, 51.6809, 51.7573, 51.8133, 51.8532, 51.9014, 51.9603,
#              52.0328, 52.0985],
#             [52.1668, 52.2316, 52.2938, 52.3680, 52.4465, 52.5180, 52.5752, 52.6178, 52.6668, 52.7340,
#              52.8056, 52.8792],
#             [52.9565, 53.0445, 53.1268, 53.2197, 53.3024, 53.3747, 53.4335, 53.4778, 53.5300, 53.5845,
#              53.6523, 53.7256],
#             [53.7882, 53.8367, 53.8830, 53.9443, 54.0042, 54.0536, 54.0856, 54.1084, 54.1463, 54.1914,
#              54.2452, 54.2958],
#             [54.3427, 54.3911, 54.4320, 54.4898, 54.5456, 54.5977, 54.6355, 54.6532, 54.6776, 54.7174,
#              54.7741, 54.8253],
#             [54.8713, 54.9161, 54.9581, 54.9997, 55.0476, 55.0912, 55.1132, 55.1328, 55.1532, 55.1898,
#              55.2416, 55.2838],
#             [55.3222, 55.3613, 55.4063, 55.4629, 55.5111, 55.5524, 55.5812, 55.6004, 55.6262, 55.6656,
#              55.7168, 55.7698],
#             [55.8197, 55.8615, 55.9130, 55.9663, 56.0220, 56.0700, 56.0939, 56.1105, 56.1314, 56.1611,
#              56.2068, 56.2583],
#             [56.3000, 56.3399, 56.3790, 56.4283, 56.4804, 56.5352, 56.5697, 56.5983, 56.6328, 56.6739,
#              56.7332, 56.7972],
#             [56.8553, 56.9111, 56.9755, 57.0471, 57.1136, 57.1738, 57.2226, 57.2597, 57.3073, 57.3643,
#              57.4334, 57.5016],
#             [57.5653, 57.6333, 57.6973, 57.7711, 57.8407, 57.9058, 57.9576, 57.9975, 58.0426, 58.1043,
#              58.1679, 58.2389],
#             [58.3092, 58.3833, 58.4537, 58.5401, 58.6228, 58.6917, 58.7410, 58.7836, 58.8406, 58.8986,
#              58.9714, 59.0438],
#             [59.1218, 59.2003, 59.2747, 59.3574, 59.4434, 59.5242, 59.5850, 59.6344, 59.6928, 59.7588,
#              59.8386, 59.9111],
#             [59.9845, 60.0564, 60.1231, 60.2042, 60.2804, 60.3530, 60.4012, 60.4440, 60.4900, 60.5578,
#              60.6324, 60.7059],
#             [60.7853, 60.8664, 60.9387, 61.0277, 61.1103, 61.1870, 61.2454, 61.2881, 61.3378, 61.4036,
#              61.4760, 61.5525],
#             [61.6287, 61.6846, 61.7433, 61.8132, 61.8823, 61.9497, 61.9969, 62.0343, 62.0714, 62.1202,
#              62.1810, 62.2382],
#             [62.2950, 62.3506, 62.3995, 62.4754, 62.5463, 62.6136, 62.6571, 62.6942, 62.7383, 62.7926,
#              62.8567, 62.9146],
#             [62.9659, 63.0217, 63.0807, 63.1462, 63.2053, 63.2599, 63.2844, 63.2961, 63.3126, 63.3422,
#              63.3871, 63.4339],
#             [63.4673, 63.4979, 63.5319, 63.5679, 63.6104, 63.6444, 63.6642, 63.6739, 63.6926, 63.7147,
#              63.7518, 63.7927],
#             [63.8285, 63.8557, 63.8804, 63.9075, 63.9393, 63.9691, 63.9799, 63.9833, 63.9938, 64.0093,
#              64.0400, 64.0670],
#             [64.0908, 64.1068, 64.1282, 64.1584, 64.1833, 64.2094, 64.2117, 64.2073, 64.2116, 64.2223,
#              64.2500, 64.2761],
#             [64.2998, 64.3192, 64.3450, 64.3735, 64.3943, 64.4151, 64.4132, 64.4118, 64.4097, 64.4168,
#              64.4329, 64.4511],
#             [64.4734, 64.4893, 64.5053, 64.5269, 64.5471, 64.5597, 64.5512, 64.5371, 64.5359, 64.5415,
#              64.5544, 64.5654],
#             [64.5736, 64.5891, 64.6015, 64.6176, 64.6374, 64.6549, 64.6530, 64.6379, 64.6372, 64.6400,
#              64.6543, 64.6723],
#             [64.6876, 64.7052, 64.7313, 64.7575, 64.7811, 64.8001, 64.7995, 64.7876, 64.7831, 64.7921,
#              64.8096, 64.8311],
#             [64.8452, 64.8597, 64.8850, 64.9175, 64.9480, 64.9794, 64.9895, 65.0028, 65.0138, 65.0371,
#              65.0773, 65.1122],
#             [65.1464, 65.1833, 65.2145, 65.2494, 65.2921, 65.3279, 65.3413, 65.3452, 65.3496, 65.3711,
#              65.3972, 65.4296],
#             [65.4573, 65.4868, 65.5152, 65.5450, 65.5781, 65.6127, 65.6288, 65.6370, 65.6493, 65.6760,
#              65.7097, 65.7461],
#             [65.7768, 65.8025, 65.8237, 65.8595, 65.8973, 65.9323, 65.9509, 65.9534, 65.9628, 65.9839,
#              66.0147, 66.0420],
#             [66.0699, 66.0961, 66.1310, 66.1683, 66.2072, 66.2356, 66.2409, 66.2335, 66.2349, 66.2441,
#              66.2751, 66.3054],
#             [66.3246, 66.3406, 66.3624, 66.3957, 66.4289, 66.4619, 66.4749, 66.4751, 66.4829, 66.5056,
#              66.5383, 66.5706],
#             [66.6030, 66.6340, 66.6569, 66.6925, 66.7289, 66.7579, 66.7708, 66.7740, 66.7846, 66.8103,
#              66.8400, 66.8779],
#             [66.9069, 66.9443, 66.9763, 67.0258, 67.0716, 67.1100, 67.1266, 67.1331, 67.1458, 67.1718,
#              67.2091, 67.2460],
#             [67.2810, 67.3136, 67.3457, 67.3890]]
#
#     def get_leap_seconds(self):
#         """returns adjustment to be added to UTC at the specified datetime to produce TAI.
#
#         Warning: This doesn't know about leap seconds after 2014 - cos I haven't worked them out ... yet
#         """
#
#         adj = 10  # as decreed from 1972
#         year = self.leap_seconds_base_year
#         while True :
#             entry = self.leap_seconds_adjustments[year - self.leap_seconds_base_year]
#             if year == self.dtt.tm_year :
#                 if self.dtt.tm_mon > 6 :
#                     adj += entry[0]
#                 break
#             adj += entry[0] + entry[1]
#             year += 1
#         return adj
#
#     def get_delta_t(self):
#         """returns a suitable value for delta_t for the given datetime."""
#
#         year, month = self.dtt.tm_year, self.dtt.tm_mon
#         if year < self.delta_t_base_year:
#             year = self.delta_t_base_year
#             month = 1
#         elif year == self.delta_t_base_year:
#             month = max(0, month - self.delta_t_base_month) + 1
#         elif year >= self.delta_t_base_year + len(self.delta_t):
#             year = self.delta_t_base_year + len(self.delta_t) - 1
#         if year == self.delta_t_base_year + len(self.delta_t) - 1:
#             month = min(month, len(self.delta_t[year - self.delta_t_base_year]))
#         return self.delta_t[year - self.delta_t_base_year][month - 1]
#
#     def get_julian_solar_day(self):
#         """returns the UT Julian day number (including fraction of a day) corresponding to the specified date/time.
#
#         This version assumes the proleptic Gregorian calendar; trying to adjust for pre-Gregorian dates/times seems
#         pointless when the changeover happened over such wildly varying times in different regions.
#         """
#         return ((self.dt.timestamp() + self.get_leap_seconds() + self.tt_offset - self.get_delta_t()) / self.seconds_per_day + self.gregorian_day_offset + self.julian_day_offset)
#
#     def get_julian_ephemeris_day(self):
#         """returns the TT Julian day number (including fraction of a day) corresponding to the specified date/time.
#
#         This version assumes the proleptic Gregorian calendar; trying to adjust for pre-Gregorian dates/times seems
#         pointless when the changeover happened over such wildly varying times in different regions.
#         """
#         return ((self.dt.timestamp() + self.get_leap_seconds() + self.tt_offset) / self.seconds_per_day + self.gregorian_day_offset + self.julian_day_offset)
#
#     def get_julian_century(self):
#         return (self.get_julian_solar_day() - 2451545.0) / 36525.0
#
#     def get_julian_ephemeris_century(self):
#         return (self.get_julian_ephemeris_day() - 2451545.0) / 36525.0
#
#     def get_julian_ephemeris_millennium(self):
#         return (self.get_julian_ephemeris_century() / 10.0)
#
#     def equation_of_time(self):
#         """returns the number of minutes to add to mean solar time to get actual solar time."""
#
#         b = 2 * self.math.pi / 364.0 * (self.yday - 81)
#         return 9.87 * self.math.sin(2 * b) - 7.53 * self.math.cos(b) - 1.5 * self.math.sin(b)
#
#     def get_coefficients(self, coeffs):
#         """computes a polynomial with time-varying coefficients from the given constant coefficients array and the current Julian millennium."""
#
#         result = 0.0
#         x = 1.0
#         for line in coeffs:
#             c = 0.0
#             for l in line:
#                 c += l[0] * self.math.cos(l[1] + l[2] * self.get_julian_ephemeris_millennium())
#             result += c * x
#             x *= self.get_julian_ephemeris_millennium()
#         return result
#
#     def get_sun_earth_distance(self):
#         return self.get_coefficients(self.sun_earth_distance_coefficients) / 1e8
#
#     def get_aberration_correction(self):
#         """sun-earth distance is in astronomical units"""
#
#         return -20.4898/(3600.0 * self.get_sun_earth_distance())
#
#     def get_equatorial_horizontal_parallax(self):
#         return 8.794 / (3600 / self.get_sun_earth_distance())
#
#     def get_nutation(self):
#         jce = self.get_julian_ephemeris_century()
#         abcd = self.nutation_coefficients
#         nutation_long = []
#         nutation_oblique = []
#         p = self.aberration_coefficients
#         x = list(p[k](jce) for k in (
#                     'MeanElongationOfMoon',
#                     'MeanAnomalyOfSun',
#                     'MeanAnomalyOfMoon',
#                     'ArgumentOfLatitudeOfMoon',
#                     'LongitudeOfAscendingNode',
#         ))
#         y = self.aberration_sin_terms
#         for i in range(len(abcd)):
#             sigmaxy = 0.0
#             for j in range(len(x)):
#                 sigmaxy += x[j] * y[i][j]
#             nutation_long.append((abcd[i][0] + (abcd[i][1] * jce)) * self.math.sin(self.math.radians(sigmaxy)))
#             nutation_oblique.append((abcd[i][2] + (abcd[i][3] * jce)) * self.math.cos(self.math.radians(sigmaxy)))
#
#         # 36000000 scales from 0.0001 arcseconds to degrees
#         nutation = {'longitude': sum(nutation_long)/36000000.0, 'obliquity': sum(nutation_oblique)/36000000.0}
#
#         return nutation
#
