from climate.constants import *
from climate.helpers import *

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from matplotlib import cm

from psychrolib import GetHumRatioFromRelHum
from windrose import WindroseAxes

import pandas as pd
import numpy as np


def plot_heatmap_new(pandas_series, cmap='Greys', tone_color="k", title=None, vrange=None):

    # Instantiate figure
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))

    # Load data
    pandas_series = pandas_series.to_frame()

    # Remove timezone-awareness
    try:
        pandas_series.index = pandas_series.index.tz_convert(None)
    except Exception as e:
        raise e

    # Reshape data into time/day matrix
    ll = pandas_series.pivot_table(columns=pandas_series.index.date, index=pandas_series.index.time).values[::-1]

    # Plot data
    heatmap = ax.imshow(
        ll,
        extent=[mdates.date2num(pandas_series.index.min()), mdates.date2num(pandas_series.index.max()), 726449, 726450],
        aspect='auto',
        cmap=cmap,
        interpolation='none',
        vmin=vrange[0] if vrange is not None else None,
        vmax=vrange[-1] if vrange is not None else None,
    )

    # Axis formatting
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.yaxis_date()
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.invert_yaxis()
    ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(ax.get_xticklabels(), ha='left', color=tone_color)
    plt.setp(ax.get_yticklabels(), color=tone_color)

    # Spine formatting
    [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]

    # Grid formatting
    ax.grid(b=True, which='major', color='white', linestyle=':', alpha=1)

    # Colorbar formatting
    cb = fig.colorbar(heatmap, orientation='horizontal', drawedges=False, fraction=0.05, aspect=100, pad=0.075)
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=tone_color)
    cb.outline.set_visible(False)

    # Add title if provided
    if title is not None:
        plt.title(title, color=tone_color, y=1.01)

    # Tidy plot
    plt.tight_layout()

    # # Save figure
    # if save:
    #     save_path.parent.mkdir(parents=True, exist_ok=True)
    #     fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
    #     print("Heatmap saved to {}".format(save_path))
    # if close:
    #     plt.close()

    return fig

def plot_heatmap(weather_object, variable, cmap='Greys', close=True, save=False):
    # Construct the save_path and create directory if it doesn't exist
    save_path = weather_object.file_path.parent / "{}_Plot".format(weather_object.file_path.stem) / "heatmap_{}.png".format(variable)

    # Instantiate figure
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))

    # Load data
    series = weather_object.df[variable].to_frame()

    # Remove timezone-awareness
    try:
        series.index = series.index.tz_convert(None)
    except Exception as e:
        raise e

    # Reshape data into time/day matrix
    ll = series.pivot_table(columns=series.index.date, index=series.index.time).values[::-1]

    # Plot data
    heatmap = ax.imshow(
        ll,
        extent=[mdates.date2num(series.index.min()), mdates.date2num(series.index.max()), 726449, 726450],
        aspect='auto',
        cmap=cmap,
        interpolation='none'
    )

    # Formatting
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.yaxis_date()
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.invert_yaxis()
    ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(ax.get_xticklabels(), ha='left', color='k')
    plt.setp(ax.get_yticklabels(), color='k')
    [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]
    ax.grid(b=True, which='major', color='white', linestyle='-', alpha=1)
    cb = fig.colorbar(heatmap, orientation='horizontal', drawedges=False, fraction=0.05, aspect=100, pad=0.075)
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='k')
    cb.outline.set_visible(False)
    plt.title("{}\n{} - {} - {}".format(renamer(series.columns[0]), weather_object.city, weather_object.country, weather_object.station_id),
              color='k', y=1.01)
    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("Heatmap saved to {}".format(save_path))
    if close:
        plt.close()

    return fig

def plot_diurnal(weather_object, dew_point=False, close=True, save=False):
    # Construct the save_path and create directory if it doesn't exist
    save_path = weather_object.file_path.parent / "{}_Plot".format(weather_object.file_path.stem) / "diurnal_{}.png".format(
        "dpt" if dew_point else "rh")

    # Group dry-bulb temperatures
    a_gp = weather_object.dry_bulb_temperature.groupby([weather_object.index.month, weather_object.index.hour])
    a_min = a_gp.min().reset_index(drop=True)
    a_mean = a_gp.mean().reset_index(drop=True)
    a_max = a_gp.max().reset_index(drop=True)

    # Group relative humidity / dewpoint temperature
    if dew_point:
        b_gp = weather_object.dew_point_temperature.groupby([weather_object.index.month, weather_object.index.hour])
        b_var = "Dew-point Temperature (°C)"
    else:
        b_gp = weather_object.relative_humidity.groupby(
            [weather_object.relative_humidity.index.month, weather_object.relative_humidity.index.hour])
        b_var = 'Relative Humidity (%)'
    b_min = b_gp.min().reset_index(drop=True)
    b_mean = b_gp.mean().reset_index(drop=True)
    b_max = b_gp.max().reset_index(drop=True)

    # Group solar radiation
    c_global_mean = weather_object.global_horizontal_radiation.groupby(
        [weather_object.index.month, weather_object.index.hour]).mean().reset_index(drop=True)
    c_diffuse_mean = weather_object.diffuse_horizontal_radiation.groupby(
        [weather_object.index.month, weather_object.index.hour]).mean().reset_index(drop=True)
    c_direct_mean = weather_object.direct_normal_radiation.groupby([weather_object.index.month, weather_object.index.hour]).mean().reset_index(
        drop=True)

    # Instantiate plot
    fig, ax = plt.subplots(3, 1, figsize=(15, 8))

    # Plot DBT
    [ax[0].plot(a_mean.iloc[i:i + 24], color='#BC204B', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[0].fill_between(np.arange(i, i + 24), a_min.iloc[i:i + 24], a_max.iloc[i:i + 24], color='#BC204B',
                        alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
    ax[0].set_ylabel('Dry-bulb Temperature (°C)', labelpad=2, color='k')
    ax[0].yaxis.set_major_locator(MaxNLocator(7))

    # Plot DPT / RH
    [ax[1].plot(b_mean.iloc[i:i + 24], color='#00617F', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[1].fill_between(np.arange(i, i + 24), b_min.iloc[i:i + 24], b_max.iloc[i:i + 24], color='#00617F',
                        alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
    ax[1].set_ylabel(b_var, labelpad=2, color='k')
    ax[1].yaxis.set_major_locator(MaxNLocator(7))
    if not dew_point:
        ax[1].set_ylim([0, 100])

    # Plot solar
    [ax[2].plot(c_direct_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls='--', label='Direct Normal Radiation') for
     i in
     np.arange(0, 288)[::24]]
    [ax[2].plot(c_diffuse_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls=':',
                label='Diffuse Horizontal Radiation') for i
     in np.arange(0, 288)[::24]]
    [ax[2].plot(c_global_mean.iloc[i:i + 24], color='#FF8F1C', lw=2, ls='-', label='Global Horizontal Radiation')
     for
     i in np.arange(0, 288)[::24]]
    ax[2].set_ylabel('Solar Radiation (W/m$^{2}$)', labelpad=2, color='k')
    ax[2].yaxis.set_major_locator(MaxNLocator(7))

    # Format plot area
    [[i.spines[spine].set_visible(False) for spine in ['top', 'right']] for i in ax]
    [[i.spines[j].set_color('k') for i in ax] for j in ['bottom', 'left']]
    [i.xaxis.set_ticks(np.arange(0, 288, 24)) for i in ax]
    [i.set_xlim([0, 287]) for i in ax]
    [plt.setp(i.get_yticklabels(), color='k') for i in ax]
    [i.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
                       ha='left', color='k') for i in ax]
    [i.get_xaxis().set_ticklabels([]) for i in [ax[0], ax[1]]]
    [i.grid(b=True, which='major', axis='both', c='k', ls='--', lw=1, alpha=0.3) for i in ax]
    [i.tick_params(length=0) for i in ax]

    ax[2].set_ylim([0, ax[2].get_ylim()[1]])

    # Legend
    handles, labels = ax[2].get_legend_handles_labels()
    lgd = ax[2].legend(bbox_to_anchor=(0.5, -0.2), loc=8, ncol=3, borderaxespad=0., frameon=False,
                       handles=[handles[0], handles[12], handles[24]], labels=[labels[0], labels[12], labels[24]])
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color='k') for text in lgd.get_texts()]

    # Add a title
    title = plt.suptitle(
        "Monthly average diurnal profile\n{0:} - {1:} - {2:}".format(weather_object.city, weather_object.country, weather_object.station_id),
        color='k', y=1.025)

    # Tidy plot
    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("Diurnal plot saved to {}".format(save_path))
    if close:
        plt.close()

    return fig

def plot_radiation_rose(self, save=False, close=False):

    rose_angles = np.radians(np.arange(0, 360, 10))
    radiation_rose_vectors = np.stack([-np.sin(rose_angles), np.cos(rose_angles), np.zeros(len(rose_angles))]).T

    direct_values = self.direct_sky_matrix.sum(axis=0)
    diffuse_values = self.diffuse_sky_matrix.sum(axis=0)

    direct_sky_radiation_rose_values = []
    diffuse_sky_radiation_rose_values = []
    for vec in radiation_rose_vectors:
        direct_radiation = 0
        diffuse_radiation = 0
        for patch_number, patch_vector in enumerate(self.patch_vectors):
            vector_angle = angle_between(patch_vector, vec, degrees=True)
            if vector_angle < 90:
                direct_radiation += direct_values[patch_number] * np.cos(np.radians(vector_angle))
                diffuse_radiation += diffuse_values[patch_number] * np.cos(np.radians(vector_angle))
        direct_sky_radiation_rose_values.append(direct_radiation)
        diffuse_sky_radiation_rose_values.append(diffuse_radiation)
    direct_sky_radiation_rose_values = np.array(direct_sky_radiation_rose_values)
    diffuse_sky_radiation_rose_values = np.array(diffuse_sky_radiation_rose_values)
    total_sky_radiation_rose_values = direct_sky_radiation_rose_values + diffuse_sky_radiation_rose_values

    max_val = total_sky_radiation_rose_values.max() / 1000
    figs = []
    for tx, vals in {"Direct Radiation": direct_sky_radiation_rose_values,
                     "Diffuse Radiation": diffuse_sky_radiation_rose_values,
                     "Total Radiation": total_sky_radiation_rose_values}.items():

        # Construct the save_path and create directory if it doesn't exist
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "radiationrose_{}.png".format(
            tx)

        vals /= 1000
        colors = [cm.OrRd(i) for i in np.interp(vals, [min(vals), max(vals)], [0, 1])]
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), subplot_kw={'projection': "polar"})
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.bar(
            rose_angles[::-1],
            vals,
            width=(np.pi * 2 / 37), zorder=5, bottom=0.0, color=colors, alpha=1, edgecolor="w",
            linewidth=0)  # 37 to make the width slightly smaller
        ax.set_ylim(0, max_val)
        ax.spines['polar'].set_visible(False)
        ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
        ti = ax.set_title("{}\n{} - {} - {}".format(tx, self.city, self.country, self.station_id), color="k",
                          loc="left", va="bottom", ha="left", fontsize="large", y=1)
        plt.tight_layout()

        # Save figure
        if save:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
            print("Radiation rose saved to {}".format(save_path))
        if close:
            plt.close()
        figs.append(fig)

    return figs

def plot_windrose(self, season_period="Annual", day_period="Daily", nsector=16, cmap=None, save=False, close=False):

    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "windrose_{}_{}.png".format(
        season_period, day_period)

    # Descibe a set of masks to remove unwanted hours of the year
    toy_masks = {
        "Daily": ((self.index.hour >= 0) & (self.index.hour <= 24)),
        "Morning": ((self.index.hour >= 5) & (self.index.hour <= 10)),
        "Midday": ((self.index.hour >= 11) & (self.index.hour <= 13)),
        "Afternoon": ((self.index.hour >= 14) & (self.index.hour <= 18)),
        "Evening": ((self.index.hour >= 19) & (self.index.hour <= 22)),
        "Night": ((self.index.hour >= 23) | (self.index.hour <= 4)),

        "Annual": ((self.index.month >= 1) & (self.index.month <= 12)),
        "Spring": ((self.index.month >= 3) & (self.index.month <= 5)),
        "Summer": ((self.index.month >= 6) & (self.index.month <= 8)),
        "Autumn": ((self.index.month >= 9) & (self.index.month <= 11)),
        "Winter": ((self.index.month <= 2) | (self.index.month >= 12))
    }
    speed_mask = (self.wind_speed != 0)
    direction_mask = (self.wind_direction != 0)
    mask = np.array([toy_masks[day_period], toy_masks[season_period], speed_mask, direction_mask]).all(axis=0)

    fig = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    unit = "m/s"
    ax.bar(self.wind_direction[mask], self.wind_speed[mask], normed=True,
           bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], opening=1, edgecolor='White',
           lw=0.25, nsector=nsector,
           cmap=plt.cm.Purples if cmap is None else cmap)

    lgd = ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', frameon=False, title=unit)
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color='k') for text in lgd.get_texts()]
    plt.setp(lgd.get_title(), color='k')

    for i, leg in enumerate(lgd.get_texts()):
        b = leg.get_text().replace('[', '').replace(')', '').split(' : ')
        lgd.get_texts()[i].set_text(b[0] + ' to ' + b[1])

    ax.grid(linestyle=':', color='k', alpha=0.5)
    ax.spines['polar'].set_visible(False)
    plt.setp(ax.get_xticklabels(), color='k')
    plt.setp(ax.get_yticklabels(), color='k')
    ax.set_title("{2:} - {3:} - {4:}\n{0:} - {1:} - {5:}".format(self.city, self.country, season_period, day_period,
                                                                 "Wind speed", self.station_id), y=1.06, color='k')

    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300, transparent=False)
        print("Windrose saved to {}".format(save_path))
    if close:
        plt.close()

    return fig

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
