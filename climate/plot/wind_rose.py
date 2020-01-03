from climate.common.constants import time_of_year_mask

import matplotlib.pyplot as plt
from matplotlib import cm
from windrose import WindroseAxes
import numpy as np


def windrose(self, season_period="Annual", day_period="Daily", n_sector=16, cmap=None, tone_color="k", save=False, close=False):

    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "windrose_{}{}.png".format(
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
           cmap=cm.get_cmap('Purples') if cmap is None else cm.get_cmap(cmap))

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
    if close:
        plt.close()

    return fig
