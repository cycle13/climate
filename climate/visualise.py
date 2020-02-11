
import numpy as np
import pandas as pd

from .weather import Weather
import pathlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from matplotlib import cm
from windrose import WindroseAxes
import colour

#############
# CONSTANTS #
#############

class Visuals(object):

    def __init__(self):

        # Colours
        self.c_bh_lightred = colour.rgb2hex([188, 32, 75], force_long=False)
        self.c_bh_brightred = colour.rgb2hex([213, 0, 50], force_long=False)
        self.c_bh_lightblue = colour.rgb2hex([141, 185, 202], force_long=False)
        self.c_bh_brightgreen = colour.rgb2hex([108, 194, 78], force_long=False)
        self.c_bh_teal = colour.rgb2hex([0, 164, 153], force_long=False)
        self.c_bh_orange = colour.rgb2hex([255, 143, 28], force_long=False)
        self.c_bh_darkblue = colour.rgb2hex([0, 97, 127], force_long=False)
        self.c_bh_green = colour.rgb2hex([34, 136, 72], force_long=False)
        self.c_bh_purple = colour.rgb2hex([36, 19, 95], force_long=False)
        self.c_bh_darkred = colour.rgb2hex([126, 45, 64], force_long=False)
        self.c_bh_grey = colour.rgb2hex([150, 140, 131], force_long=False)
        self.c_bh_brightpurple = colour.rgb2hex([112, 47, 138], force_long=False)
        self.c_bh_pink = colour.rgb2hex([231, 130, 169], force_long=False)
        self.c_bh_brightblue = colour.rgb2hex([0, 169, 224], force_long=False)
        self.c_bh_purple_ppt = colour.rgb2hex([143, 114, 176], force_long=False)
        self.c_bh_brandgreen = colour.rgb2hex([196, 214, 0], force_long=False)
        self.c_bh_pink_ppt = colour.rgb2hex([230, 49, 135], force_long=False)
        self.c_bh_blue = colour.rgb2hex([0, 176, 240], force_long=False)
        self.c_bh_darkblue_ppt = colour.rgb2hex([28, 54, 96], force_long=False)
        self.c_bh_orange_ppt = colour.rgb2hex([235, 103, 28], force_long=False)
        self.c_bh_darkgrey = colour.rgb2hex([63, 63, 63], force_long=False)

        # Heatmaps
        self.cm_wind_speed = cm.get_cmap("GnBu")
        self.cm_wind_direction = cm.get_cmap("PuBu")
        self.cm_direct_normal_radiation = cm.get_cmap("OrRd")
        self.cm_diffuse_horizontal_radiation = cm.get_cmap("OrRd")
        self.cm_global_horizontal_radiation = cm.get_cmap("OrRd")
        self.cm_dry_bulb_temperature = cm.get_cmap("Reds")
        self.cm_relative_humidity = cm.get_cmap("Blues")
        self.cm_dew_point_temperature = cm.get_cmap("Greens")

        utci = ListedColormap([item for sublist in [
            ['#1A5899'] * 13,
            ['#347FB9'] * 14,
            ['#82BBD9'] * 13,
            ['#BFDCEB'] * 10,
            ['#FFFFFF'] * 17,
            ['#F7C1AA'] * 6,
            ['#E3806B'] * 6,
            ['#C84648'] * 8
        ] for item in sublist])
        utci.set_under('#053061')
        utci.set_over('#B2182B')
        self.cm_utci = utci

        utci_ghadan = ListedColormap([item for sublist in [
            ['#262972'] * 13,
            ['#3452A4'] * 14,
            ['#3C65AF'] * 13,
            ['#37BCED'] * 10,
            ['#2EB349'] * 17,
            ['#F38322'] * 6,
            ['#C31F25'] * 6,
            ['#7F1416'] * 8
        ] for item in sublist])
        utci_ghadan.set_under('#0D104B')
        utci_ghadan.set_over('#580002')
        self.cm_utci_ghadan = utci_ghadan


class TimePeriod(object):

    def __init__(self, DATETIME_INDEX = pd.date_range(start="2018-01-01 00:00:00", end="2019-01-01 00:00:00", freq="60T", closed="left")):

        self.daily = ((DATETIME_INDEX.hour >= 0) & (DATETIME_INDEX.hour <= 24))
        self.morning = ((DATETIME_INDEX.hour >= 5) & (DATETIME_INDEX.hour <= 10))
        self.midday = ((DATETIME_INDEX.hour >= 11) & (DATETIME_INDEX.hour <= 13))
        self.afternoon = ((DATETIME_INDEX.hour >= 14) & (DATETIME_INDEX.hour <= 18))
        self.evening = ((DATETIME_INDEX.hour >= 19) & (DATETIME_INDEX.hour <= 22))
        self.night = ((DATETIME_INDEX.hour >= 23) | (DATETIME_INDEX.hour <= 4))
        self.shoulder_morning = ((DATETIME_INDEX.hour >= 7) & (DATETIME_INDEX.hour <= 10))
        self.shoulder_afternoon = ((DATETIME_INDEX.hour >= 16) & (DATETIME_INDEX.hour <= 19))
        self.annual = ((DATETIME_INDEX.month >= 1) & (DATETIME_INDEX.month <= 12))
        self.spring = ((DATETIME_INDEX.month >= 3) & (DATETIME_INDEX.month <= 5))
        self.summer = ((DATETIME_INDEX.month >= 6) & (DATETIME_INDEX.month <= 8))
        self.autumn = ((DATETIME_INDEX.month >= 9) & (DATETIME_INDEX.month <= 11))
        self.winter = ((DATETIME_INDEX.month <= 2) | (DATETIME_INDEX.month >= 12))
        self.shoulder_month = ((DATETIME_INDEX.month == 3) | (DATETIME_INDEX.month == 10))


def wind_rose(weather, day_period="daily", season_period="annual", n_sector=16, cmap=Visuals().cm_wind_speed, tone_color="#555555", save=False, close=True):

    # Describe a set of masks to remove unwanted hours of the year
    mask = np.array([getattr(TimePeriod(), day_period), getattr(TimePeriod(), season_period), (weather.wind_speed != 0), (weather.wind_direction != 0)]).all(axis=0)

    fig = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    ax.bar(
        weather.wind_direction[mask], weather.wind_speed[mask],
        normed=True,
        bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        opening=0.95,
        edgecolor='White',
        lw=0.1,
        nsector=n_sector,
        cmap=cmap
    )

    lgd = ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', frameon=False, title="m/s")
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color=tone_color) for text in lgd.get_texts()]
    plt.setp(lgd.get_title(), color=tone_color)

    for i, leg in enumerate(lgd.get_texts()):
        b = leg.get_text().replace('[', '').replace(')', '').split(' : ')
        lgd.get_texts()[i].set_text(b[0] + ' to ' + b[1])

    ax.grid(linestyle=':', color=tone_color, alpha=0.5)
    ax.spines['polar'].set_visible(False)
    plt.setp(ax.get_xticklabels(), color=tone_color)
    plt.setp(ax.get_yticklabels(), color=tone_color)
    ax.set_title("{} - {}\n{} - {} - {}".format(season_period.title(), day_period.title(), weather.city, weather.country, weather.station_id),
                 y=1.06, color=tone_color, loc="center", va="bottom", ha="center", fontsize="medium")

    # plt.tight_layout()

    if save:
        save_path = pathlib.Path(weather.epw_path).parent / "{}_plot".format(pathlib.Path(weather.epw_path).stem) / "wind_rose_{}{}.png".format(season_period, day_period)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_extra_artists=(lgd, ), bbox_inches='tight', dpi=300, transparent=False)
        print("Windrose saved to {}".format(save_path))

    if close:
        plt.close()

    return fig


