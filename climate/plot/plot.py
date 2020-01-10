from climate.common.helpers import *

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from matplotlib import cm

from psychrolib import GetHumRatioFromRelHum
from windrose import WindroseAxes

import numpy as np





def plot_psychrometrics(self, nbins=50, cm="Greys", close=False, save=False):
    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "psychrometrics.png"

    hr = self.humidity_ratio
    dbt = self.dry_bulb_temperature

    dry_bulb_temperature_plot_range = range(-20, 51, 1)
    humidity_ratio_plot_range = [i / 10000 for i in range(0, 301, 1)]
    enthalpy_plot_range = range(-10, 120, 10)
    relative_humidity_plot_range = [i / 100 for i in range(0, 101, 10)]

    # figure instantiation
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))

    # plot values from weather file
    counts, xedges, yedges, im = ax.hist2d(dbt, hr, bins=nbins, cmin=1, alpha=0.9, normed=False, cmap=cm, lw=0,
                                           zorder=0)

    # y-axis formatting
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylim(0, 0.03)
    ax.set_yticks([i / 1000 for i in range(0, 35, 5)])
    ax.set_ylabel("Humidity ratio ($kg_{water}/kg_{air}$)", color="k", fontsize="x-large")

    # x-axis formatting
    ax.set_xlim(-20, 50)
    ax.set_xticks(range(-20, 55, 5))
    ax.set_xlabel("Dry-bulb temperature ($°C$)", color="k", fontsize="x-large")
    ax.tick_params(axis='both', colors='k')

    # canvas formatting
    ax.tick_params(axis="both", color="k", grid_color="k", grid_alpha=1, grid_lw=0.5)
    for edge in ["right", "bottom"]:
        ax.spines[edge].set_alpha(1)
        ax.spines[edge].set_color("k")
        ax.spines[edge].set_lw(1)
    for edge in ["left", "top"]:
        ax.spines[edge].set_visible(False)

    # relative humidity grid/curves
    n = 0
    for rh in relative_humidity_plot_range:
        h_r = [GetHumRatioFromRelHum(i, rh, 101350) for i in dry_bulb_temperature_plot_range]
        ax.plot(dry_bulb_temperature_plot_range, h_r, color="k", alpha=1, lw=0.2)
        # Fill the top part of the plot
        if rh == 1:
            ax.fill_between(dry_bulb_temperature_plot_range, h_r, 0.031, interpolate=True, color='w', lw=0,
                            edgecolor=None,
                            zorder=4)
        # add annotation describing line
        ax.text(30, GetHumRatioFromRelHum(30, rh, 101350) + 0.0, "{0:0.0f}% RH".format(rh * 100),
                ha="right", va="bottom", rotation=0, zorder=9, fontsize="small", color="k")  # n*55
        n += 1 / len(relative_humidity_plot_range)

    # TODO: Fix enthalpy grid curves
    # # enthalpy grid/curves
    # for enthalpy in enthalpy_plot_range:
    #     ys = [0, 0.030]
    #     xs = np.array([GetTDryBulbFromEnthalpyAndHumRatio(enthalpy, i) for i in ys]) /1000
    #     if (enthalpy <= 50) & (enthalpy != 30):
    #         ax.text(xs[0], 0.0002, "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="k", zorder=9,
    #                 fontsize="small")
    #     else:
    #         pass
    #     # ax.text(50, ys[0], "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="#555555", zorder=9, fontsize="small")
    #     ax.plot(xs, ys, color="k", alpha=1, lw=0.2)

    # grid formatting
    ax.grid(True, lw=0.2, zorder=5)

    # Generating summary metrics
    min_dbt = self.df[self.index == self.dry_bulb_temperature.idxmin()].squeeze()
    max_dbt = self.df[self.index == self.dry_bulb_temperature.idxmax()].squeeze()
    max_hr = self.df[self.index == self.humidity_ratio.idxmax()].squeeze()
    max_enthalpy = self.df[self.index == self.enthalpy.idxmax()].squeeze()

    text_fontsize = "medium"

    # Generate peak cooling summary
    max_dbt_table = "Peak cooling {0:}\n" \
                    "WS:  {1:>6.1f} m/s\n" \
                    "WD:  {2:>6.1f} deg\n" \
                    "DBT: {3:>6.1f} °C\n" \
                    "WBT: {4:>6.1f} °C\n" \
                    "RH:  {5:>6.1f} %\n" \
                    "DPT: {6:>6.1f} °C\n" \
                    "h:   {7:>6.1f} kJ/kg\n" \
                    "HR:  {8:<5.4f} kg/kg".format(
        max_dbt.name.strftime("%b %d %H:%M"),
        max_dbt.wind_speed,
        max_dbt.wind_direction,
        max_dbt.dry_bulb_temperature,
        max_dbt.wet_bulb_temperature,
        max_dbt.relative_humidity,
        max_dbt.dew_point_temperature,
        max_dbt.enthalpy / 1000,
        max_dbt.humidity_ratio)

    ax.text(0, 0.98, max_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="k", **{'fontname': 'monospace'})

    ## Generate peak heating summary
    min_dbt_table = "Peak heating {0:}\n" \
                    "WS:  {1:>6.1f} m/s\n" \
                    "WD:  {2:>6.1f} deg\n" \
                    "DBT: {3:>6.1f} °C\n" \
                    "WBT: {4:>6.1f} °C\n" \
                    "RH:  {5:>6.1f} %\n" \
                    "DPT: {6:>6.1f} °C\n" \
                    "h:   {7:>6.1f} kJ/kg\n" \
                    "HR:  {8:<5.4f} kg/kg".format(
        min_dbt.name.strftime("%b %d %H:%M"),
        min_dbt.wind_speed,
        min_dbt.wind_direction,
        min_dbt.dry_bulb_temperature,
        min_dbt.wet_bulb_temperature,
        min_dbt.relative_humidity,
        min_dbt.dew_point_temperature,
        min_dbt.enthalpy / 1000,
        min_dbt.humidity_ratio
    )
    ax.text(0, 0.72, min_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="k", **{'fontname': 'monospace'})

    ## Generate max HumidityRatio summary
    max_hr_table = "Peak humidity ratio {0:}\n" \
                   "WS:  {1:>6.1f} m/s\n" \
                   "WD:  {2:>6.1f} deg\n" \
                   "DBT: {3:>6.1f} °C\n" \
                   "WBT: {4:>6.1f} °C\n" \
                   "RH:  {5:>6.1f} %\n" \
                   "DPT: {6:>6.1f} °C\n" \
                   "h:   {7:>6.1f} kJ/kg\n" \
                   "HR:  {8:<5.4f} kg/kg".format(
        max_hr.name.strftime("%b %d %H:%M"),
        max_hr.wind_speed,
        max_hr.wind_direction,
        max_hr.dry_bulb_temperature,
        max_hr.wet_bulb_temperature,
        max_hr.relative_humidity,
        max_hr.dew_point_temperature,
        max_hr.enthalpy / 1000,
        max_hr.humidity_ratio
    )
    ax.text(0.17, 0.98, max_hr_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="k", **{'fontname': 'monospace'})

    ## Generate max enthalpy summary
    max_enthalpy_table = "Peak enthalpy ratio {0:}\n" \
                         "WS:  {1:>6.1f} m/s\n" \
                         "WD:  {2:>6.1f} deg\n" \
                         "DBT: {3:>6.1f} °C\n" \
                         "WBT: {4:>6.1f} °C\n" \
                         "RH:  {5:>6.1f} %\n" \
                         "DPT: {6:>6.1f} °C\n" \
                         "h:   {7:>6.1f} kJ/kg\n" \
                         "HR:  {8:<5.4f} kg/kg".format(
        max_enthalpy.name.strftime("%b %d %H:%M"),
        max_enthalpy.wind_speed,
        max_enthalpy.wind_direction,
        max_enthalpy.dry_bulb_temperature,
        max_enthalpy.wet_bulb_temperature,
        max_enthalpy.relative_humidity,
        max_enthalpy.dew_point_temperature,
        max_enthalpy.enthalpy / 1000,
        max_enthalpy.humidity_ratio
    )
    ax.text(0.17, 0.72, max_enthalpy_table, transform=ax.transAxes, ha="left", va="top", zorder=8,
            fontsize=text_fontsize, color="k", **{'fontname': 'monospace'})

    # Title formatting
    ti = ax.set_title("{} - {} - {}".format(self.city, self.country, self.station_id),
                      color="k", loc="left", fontsize="xx-large")

    # text
    keys = "WS: Wind speed | WD: Wind direction | DBT: Dry-bulb temperature | WBT: Wet-bulb temperature\nRH: Relative humidity | DPT: Dew-point temperature | h: Enthalpy | HR: Humidity ratio"
    te = ax.text(0.5, -0.1, keys, transform=ax.transAxes, ha="center", va="top", zorder=8, fontsize="medium",
                 color="k", **{'fontname': 'monospace'})

    # Colorbar
    cb = plt.colorbar(im, ax=ax, shrink=1, pad=0.071)
    cb.ax.set_title('Hours', color="k")
    cb.outline.set_visible(False)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='k')

    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False, bbox_extra_artists=[ti, te, ])
        print("Psychrometric plot saved to {}".format(save_path))
    if close:
        plt.close()

# TODO: plot_wind_weibull, plot_utci_frequency, plot_utci_heatmap
