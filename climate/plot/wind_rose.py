from climate.common.constants import time_of_year_mask

from climate.common.helpers import rgb_to_hex
import matplotlib.ticker as mtick
import matplotlib as mpl

import matplotlib.pyplot as plt
from matplotlib import cm
from windrose import WindroseAxes
import numpy as np


def windrose(self, season_period="Annual", day_period="Daily", n_sector=16, cmap=None, tone_color="k", save=False):

    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "wind_rose_{}{}.png".format(
        season_period, day_period)

    # Describe a set of masks to remove unwanted hours of the year
    speed_mask = (self.wind_speed != 0)
    direction_mask = (self.wind_direction != 0)
    mask = np.array([time_of_year_mask[day_period], time_of_year_mask[season_period], speed_mask, direction_mask]).all(axis=0)

    fig = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    ax.bar(self.wind_direction[mask], self.wind_speed[mask], normed=True,
           bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], opening=0.95, edgecolor='White',
           lw=0.1, nsector=n_sector,
           cmap=cm.get_cmap('GnBu') if cmap is None else cm.get_cmap(cmap))

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
    ax.set_title("{} - {}\n{} - {} - {}".format(season_period, day_period, self.city, self.country, self.station_id), y=1.06, color=tone_color, loc="center", va="bottom", ha="center", fontsize="medium")

    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300, transparent=False)
        print("Windrose saved to {}".format(save_path))

    plt.close()

    return fig


def wind_frequency(self, season_period="Annual", day_period="Daily", tone_color="k", save=False):
    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "wind_frequency_{}{}.png".format(
        season_period, day_period)

    # Describe a set of masks to remove unwanted hours of the year
    speed_mask = (self.wind_speed != 0)
    direction_mask = (self.wind_direction != 0)
    mask = np.array([time_of_year_mask[day_period], time_of_year_mask[season_period], speed_mask, direction_mask]).all(
        axis=0)

    a = self.wind_speed[mask]

    bins = np.arange(0, 16, 1)

    fig, ax = plt.subplots(1, 1, figsize=(8, 3))
    a.plot(kind="hist", density=True, bins=bins, color=rgb_to_hex([28, 54, 96]), zorder=5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(18))
    ax.set_xlabel('Wind speed (m/s)', color=tone_color)
    ax.tick_params(axis='both', colors=tone_color)
    ax.set_ylabel('Frequency', color=tone_color)
    ax.tick_params(axis='both', which='major')
    ax.grid(b=True, which='major', color=tone_color, linestyle=':', alpha=0.5, zorder=3)
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    [ax.spines[j].set_color(tone_color) for j in ['bottom', 'left']]
    ti = plt.title(
        "{2:} - {3:}\n{0:} - {1:} - {4:}".format(self.city, self.country, season_period, day_period, self.station_id),
        color=tone_color)
    ax.set_xlim([0, 16])

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
    bars = [rect for rect in ax.get_children() if isinstance(rect, mpl.patches.Rectangle)]
    for bar in bars[:-1]:
        ax.text(bar.xy[0] + bar.get_width() / 2, bar.get_height() + 0.005, "{:0.0%}".format(bar.get_height()),
                ha="center", va="bottom", color=tone_color)

    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_extra_artists=(ti,), bbox_inches='tight', dpi=300, transparent=False)
        print("Wind frequency histogram saved to {}".format(save_path))

    plt.close()

    return fig
