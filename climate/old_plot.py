from climate.old import *

import matplotlib.pyplot as plt
from matplotlib import dates, cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
from windrose import WindroseAxes
import pandas as pd
import numpy as np

def psychrometric_chart(df, nbins=50, cm="Greys", close=False, savepath=False):
    hr = df.HumidityRatio
    dbt = df.dry_bulb_temperature

    dry_bulb_temperatures = range(-20, 51, 1)
    humidity_ratios = [i / 10000 for i in range(0, 301, 1)]
    enthalpys = range(-10, 120, 10)
    relative_humiditys = [i / 100 for i in range(0, 101, 10)]

    # figure instantiation
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))

    # plot values from weatherfile
    counts, xedges, yedges, im = ax.hist2d(dbt, hr, bins=nbins, cmin=1, alpha=0.9, normed=False, cmap=cm, lw=0,
                                           zorder=0)

    # y-axis formatting
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylim(0, 0.03)
    ax.set_yticks([i / 1000 for i in range(0, 35, 5)])
    ax.set_ylabel("Humidity ratio ($kg_{water}/kg_{air}$)", color="#555555", fontsize="x-large")

    # x-axis formatting
    ax.set_xlim(-20, 50)
    ax.set_xticks(range(-20, 55, 5))
    ax.set_xlabel("Dry-bulb temperature ($°C$)", color="#555555", fontsize="x-large")
    ax.tick_params(axis='both', colors='#555555')

    # canvas formatting
    ax.tick_params(axis="both", color="#555555", grid_color="#555555", grid_alpha=1, grid_lw=0.5)
    for edge in ["right", "bottom"]:
        ax.spines[edge].set_alpha(1)
        ax.spines[edge].set_color("#555555")
        ax.spines[edge].set_lw(1)
    for edge in ["left", "top"]:
        ax.spines[edge].set_visible(False)

    # relative humidity grid/curves
    n = 0
    for rh in relative_humiditys:
        h_r = [humidity_ratio_relative_humidity(dry_bulb_temperature=i, relative_humidity=rh, pressure=101.35) for i
               in dry_bulb_temperatures]
        ax.plot(dry_bulb_temperatures, h_r, color="#555555", alpha=1, lw=0.2)
        # Fill the top part of the plot
        if rh == 1:
            ax.fill_between(dry_bulb_temperatures, h_r, 0.031, interpolate=True, color='w', lw=0, edgecolor=None,
                            zorder=4)
        # add annotation describing line
        ax.text(30, humidity_ratio_relative_humidity(dry_bulb_temperature=30, relative_humidity=rh,
                                                       pressure=101.35) + 0.0, "{0:0.0f}% RH".format(rh * 100),
                ha="right", va="bottom", rotation=0, zorder=9, fontsize="small", color="#555555")  # n*55
        n += 1 / len(relative_humiditys)

    # enthalpy grid/curves
    for enthalpy in enthalpys:
        ys = [0, 0.030]
        xs = [temperature_enthalpy_humidity_ratio(enthalpy=enthalpy, humidity_ratio=i) for i in ys]
        if (enthalpy <= 50) & (enthalpy != 30):
            ax.text(xs[0], 0.0002, "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="#555555", zorder=9,
                    fontsize="small")
        else:
            pass
        #             ax.text(50, ys[0], "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="#555555", zorder=9, fontsize="small")
        ax.plot(xs, ys, color="#555555", alpha=1, lw=0.2)

    # grid formatting
    ax.grid(True, lw=0.2, zorder=5)

    # Generating summary metrics
    min_dbt = df.sort_values(["dry_bulb_temperature"]).iloc[0, :]
    max_dbt = df.sort_values(["dry_bulb_temperature"]).iloc[-1, :]
    max_hr = df.sort_values(["HumidityRatio"]).iloc[-1, :]
    max_enthalpy = df.sort_values(["Enthalpy"]).iloc[-1, :]

    text_fontsize = "medium"

    ## Generate peak cooling summary
    #     max_dbt_table = "Peak cooling {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(max_dbt.name.strftime("%b %d %H:%M"), max_dbt.WindSpeed, max_dbt.WindDirection, max_dbt.DryBulbTemperature, max_dbt.WetBulbTemperature, max_dbt.RelativeHumidity, max_dbt.DewPointTemperature, max_dbt.Enthalpy, max_dbt.HumidityRatio)
    max_dbt_table = "Peak cooling {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(
        "", max_dbt.wind_speed, max_dbt.wind_direction, 32.1, 20.8, 35.9, 15.1, 59.8, 0.010731)
    ax.text(0, 0.98, max_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="#555555", **{'fontname': 'monospace'})

    ## Generate peak heating summary
    #     min_dbt_table = "Peak heating {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(min_dbt.name.strftime("%b %d %H:%M"), min_dbt.WindSpeed, min_dbt.WindDirection, min_dbt.DryBulbTemperature, min_dbt.WetBulbTemperature, min_dbt.RelativeHumidity, min_dbt.DewPointTemperature, min_dbt.Enthalpy, min_dbt.HumidityRatio)
    min_dbt_table = "Peak heating {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(
        "", min_dbt.wind_speed, min_dbt.wind_direction, -4.0, -4.6, 86.9, -5.6, 1.8, 0.002342)
    ax.text(0, 0.72, min_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="#555555", **{'fontname': 'monospace'})

    ## Generate max HumidityRatio summary
    max_hr_table = "Peak humidity ratio {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(
        max_hr.name.strftime("%b %d %H:%M"), max_hr.wind_speed, max_hr.wind_direction, max_hr.dry_bulb_temperature,
        max_hr.WetBulbTemperature, max_hr.relative_humidity, max_hr.dew_point_temperature, max_hr.Enthalpy,
        max_hr.HumidityRatio)
    ax.text(0.17, 0.98, max_hr_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize,
            color="#555555", **{'fontname': 'monospace'})

    ## Generate max enthalpy summary
    max_enthalpy_table = "Peak enthalpy {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(
        max_enthalpy.name.strftime("%b %d %H:%M"), max_enthalpy.wind_speed, max_enthalpy.wind_direction,
        max_enthalpy.dry_bulb_temperature, max_enthalpy.WetBulbTemperature, max_enthalpy.relative_humidity,
        max_enthalpy.dew_point_temperature, max_enthalpy.Enthalpy, max_enthalpy.HumidityRatio)
    ax.text(0.17, 0.72, max_enthalpy_table, transform=ax.transAxes, ha="left", va="top", zorder=8,
            fontsize=text_fontsize, color="#555555", **{'fontname': 'monospace'})

    # Title formatting
    ti = ax.set_title("{} - {} - {}".format(df.City.values[0], df.Country.values[0], df.StationID.values[0]),
                      color="#555555", loc="left", fontsize="xx-large")

    # text
    keys = "WS: Wind speed | WD: Wind direction | DBT: Dry-bulb temperature | WBT: Wet-bulb temperature\nRH: Relative humidity | DPT: Dew-point temperature | h: Enthalpy | HR: Humidity ratio"
    te = ax.text(0.5, -0.1, keys, transform=ax.transAxes, ha="center", va="top", zorder=8, fontsize="medium",
                 color="#555555", **{'fontname': 'monospace'})

    # Colorbar
    cb = plt.colorbar(im, ax=ax, shrink=1, pad=0.071)
    cb.ax.set_title('Hours', color="#555555")
    cb.outline.set_visible(False)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='#555555')

    plt.tight_layout()

    if savepath:
        print("Saving to {}".format(savepath))
        fig.savefig(savepath, dpi=300, transparent=False, bbox_inches="tight", bbox_extra_artists=[ti, te, ])
    if close:
        plt.close()


def monthly_diurnal_plot(df, savepath=None, close=True):
    def group(series):
        result = [series.groupby([series.index.month, series.index.hour]).min().reset_index(drop=True),
                  series.groupby([series.index.month, series.index.hour]).mean().reset_index(drop=True),
                  series.groupby([series.index.month, series.index.hour]).max().reset_index(drop=True)]
        return result

    dbtMin, dbtMean, dbtMax = group(df.dry_bulb_temperature)
    rhMin, rhMean, rhMax = group(df.relative_humidity)
    dirSolMean = df.direct_normal_radiation.groupby([df.index.month, df.index.hour]).mean().reset_index(drop=True)
    diffSolMean = df.diffuse_horizontal_radiation.groupby([df.index.month, df.index.hour]).mean().reset_index(drop=True)
    globalHorizMean = df.global_horizontal_radiation.groupby([df.index.month, df.index.hour]).mean().reset_index(
        drop=True)

    fig, ax = plt.subplots(3, 1, figsize=(15, 8))
    [ax[0].plot(dbtMean.iloc[i:i + 24], color='#BC204B', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[0].fill_between(np.arange(i, i + 24), dbtMin.iloc[i:i + 24], dbtMax.iloc[i:i + 24], color='#BC204B', alpha=0.2,
                        label='Range') for i in np.arange(0, 288)[::24]]
    ax[0].set_ylabel('Dry-bulb temperature (°C)', labelpad=2, color='#555555')
    ax[0].yaxis.set_major_locator(MaxNLocator(7))
    [ax[1].plot(rhMean.iloc[i:i + 24], color='#00617F', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[1].fill_between(np.arange(i, i + 24), rhMin.iloc[i:i + 24], rhMax.iloc[i:i + 24], color='#00617F', alpha=0.2,
                        label='Range') for i in np.arange(0, 288)[::24]]
    ax[1].set_ylabel('Relative humidity (%)', labelpad=2, color='#555555')
    [ax[2].plot(dirSolMean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls='--', label='Direct normal radiation') for i in
     np.arange(0, 288)[::24]]
    [ax[2].plot(diffSolMean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls=':', label='Diffuse horizontal radiation') for i
     in np.arange(0, 288)[::24]]
    [ax[2].plot(globalHorizMean.iloc[i:i + 24], color='#FF8F1C', lw=2, ls='-', label='Global horizontal radiation') for
     i in np.arange(0, 288)[::24]]
    ax[2].set_ylabel('Solar radiation (W/m$^{2}$)', labelpad=2, color='#555555')
    ax[2].yaxis.set_major_locator(MaxNLocator(7))

    # Format plot area
    [[i.spines[spine].set_visible(False) for spine in ['top', 'right']] for i in ax]
    [[i.spines[j].set_color('#555555') for i in ax] for j in ['bottom', 'left']]
    [i.xaxis.set_ticks(np.arange(0, 288, 24)) for i in ax]
    [i.set_xlim([0, 287]) for i in ax]
    [plt.setp(i.get_yticklabels(), color='#555555') for i in ax]
    [i.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
                       ha='left', color='#555555') for i in ax]
    [i.get_xaxis().set_ticklabels([]) for i in [ax[0], ax[1]]]
    [i.grid(b=True, which='major', axis='both', c='#555555', ls='--', lw=1, alpha=0.3) for i in ax]
    [i.tick_params(length=0) for i in ax]
    ax[1].set_ylim([0, 100])
    ax[2].set_ylim([0, ax[2].get_ylim()[1]])

    # Legend
    handles, labels = ax[2].get_legend_handles_labels()
    lgd = ax[2].legend(bbox_to_anchor=(0.5, -0.2), loc=8, ncol=3, borderaxespad=0., frameon=False,
                       handles=[handles[0], handles[12], handles[24]], labels=[labels[0], labels[12], labels[24]])
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color='#555555') for text in lgd.get_texts()]

    # Add a title
    title = plt.suptitle(
        "Monthly average diurnal profile\n{0:} - {1:} - {2:}".format(df.City[0], df.Country[0], df.StationID[0]),
        color='#555555', y=1.025)

    # Tidy plot
    plt.tight_layout()

    # Save figure
    if savepath:
        print("Saving to {}".format(savepath))
        fig.savefig(savepath, bbox_extra_artists=[title, lgd, ], bbox_inches="tight", dpi=450, transparent=False)
    if close:
        plt.close()


def weather_heatmap(df, variable, cmap='Greys', savepath=None, close=True):
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    series = df[variable].to_frame()
    try:
        series.index = series.index.tz_convert(None)
    except:
        pass

    ll = series.pivot_table(columns=series.index.date, index=series.index.time).values[::-1]
    heatmap = ax.imshow(ll,
                        extent=[dates.date2num(series.index.min()), dates.date2num(series.index.max()), 726449, 726450],
                        aspect='auto', cmap=cmap, interpolation='none')
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    ax.yaxis_date()
    ax.yaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax.invert_yaxis()
    ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(ax.get_xticklabels(), ha='left', color='#555555')
    plt.setp(ax.get_yticklabels(), color='#555555')
    [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]
    ax.grid(b=True, which='major', color='white', linestyle='-', alpha=1)
    cb = fig.colorbar(heatmap, orientation='horizontal', drawedges=False, fraction=0.05, aspect=100, pad=0.075)
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='#555555')
    cb.outline.set_visible(False)
    plt.title(
        "{}\n{} - {} - {}".format(series.columns[0], df.City.values[0], df.Country.values[0], df.StationID.values[0]),
        color='#555555', y=1.01)
    plt.tight_layout()

    # Save figure
    if savepath != None:
        print("Saving to {}".format(savepath))
        fig.savefig(savepath, bbox_inches="tight", dpi=300, transparent=False)
    if close:
        plt.close()


def utci_frequency(df, hours=['00:00', '23:59'], savepath=None, close=True):
    series = df.UniversalThermalClimateIndex
    series = series.rolling(2).mean()
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    data = series.between_time(hours[0], hours[1], include_start=True, include_end=False)
    data.plot.hist(bins=np.arange(-40, 58, 1), ax=ax, zorder=4, color='#555555', alpha=0.8, density=True)
    ax.set_yticklabels(['{:0.1f}%'.format(x * 100) for x in ax.get_yticks()])
    ax.set_xlabel('UTCI (C)', color='#555555', labelpad=20)
    ax.tick_params(axis='both', colors='#555555')
    ax.set_ylabel('Frequency', color='#555555')
    ax.set_xticks([-40, -27, -13, 0, 9, 26, 32, 38, 46])
    ax.tick_params(axis='both', which='major')
    ax.grid(b=True, which='major', color='white', linestyle='--', alpha=0.9, zorder=10)
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    [ax.spines[j].set_color('#555555') for j in ['bottom', 'left']]
    plt.title('UTCI approximation (°C) annual frequency between {2:} and {3:}\n{0:} - {1:} - {4:}'.format(df.City[0],
                                                                                                          df.Country[0],
                                                                                                          hours[0],
                                                                                                          hours[1],
                                                                                                          df.StationID.values[
                                                                                                              0]),
              color='#555555', y=1.05)
    ax.axvspan(-50, -40, color='#053061', zorder=1)
    ax.axvspan(-40, -27, color='#1A5899', zorder=1)
    ax.axvspan(-27, -13, color='#347FB9', zorder=1)
    ax.axvspan(-13, 0, color='#82BBD9', zorder=1)
    ax.axvspan(0, 9, color='#BFDCEB', zorder=1)
    ax.axvspan(9, 26, color='#FFFFFF', zorder=1)
    ax.axvspan(26, 32, color='#F7C1AA', zorder=1)
    ax.axvspan(32, 38, color='#E3806B', zorder=1)
    ax.axvspan(38, 46, color='#C84648', zorder=1)
    ax.axvspan(46, 60, color='#B2182B', zorder=1)
    ax.set_xlim([-50, 60])
    bottom, top = ax.get_ylim()
    for i, j, k in zip([-45, -33.5, -20, -6.5, 4.5, 17.5, 29, 35, 42, 53], ['{0:.1f}%'.format(i) for i in
                                                                            np.histogram(data,
                                                                                         bins=[-100, -40, -27, -13, 0,
                                                                                               9, 26, 32, 38, 46, 100])[
                                                                                0].astype("float") / len(data) * 100],
                       ["Extreme\ncold stress", "Very strong\ncold stress", "Strong\ncold stress",
                        "Moderate\ncold stress", "Slight\ncold stress", "No thermal stress", "Moderate\nheat\nstress",
                        "Strong\nheat\nstress", "Very strong\nheat stress", "Extreme\nheat stress"]):
        ax.text(i, top + 0.00025, j, ha='center', va='bottom', color='#555555')
        ax.text(i, bottom - 0.001, k, ha='center', va='top', color='#555555', fontsize='small')
    plt.tight_layout()

    # Save figure
    if savepath != None:
        print("Saving to {}".format(savepath))
        fig.savefig(savepath, bbox_inches="tight", dpi=300, transparent=False)
    if close:
        plt.close()


def utci_heatmap(df, y_move=-0.22, savepath=None, close=True):
    series = df.UniversalThermalClimateIndex.to_frame()
    series = series.rolling(2).mean()
    colors = ['#1A5899', '#1A5899', '#1A5899', '#1A5899', '#1A5899', '#1A5899', '#1A5899', '#1A5899', '#1A5899',
              '#1A5899', '#1A5899', '#1A5899', '#1A5899', '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9',
              '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9', '#347FB9',
              '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9',
              '#82BBD9', '#82BBD9', '#82BBD9', '#82BBD9', '#BFDCEB', '#BFDCEB', '#BFDCEB', '#BFDCEB', '#BFDCEB',
              '#BFDCEB', '#BFDCEB', '#BFDCEB', '#BFDCEB', '#BFDCEB', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF',
              '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF',
              '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#F7C1AA', '#F7C1AA', '#F7C1AA', '#F7C1AA', '#F7C1AA',
              '#F7C1AA', '#E3806B', '#E3806B', '#E3806B', '#E3806B', '#E3806B', '#E3806B', '#C84648', '#C84648',
              '#C84648', '#C84648', '#C84648', '#C84648', '#C84648', '#C84648']
    cmap = ListedColormap(colors)
    cmap.set_under('#053061')
    cmap.set_over('#B2182B')
    bounds = np.arange(-41, 48, 1)
    norm = BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    heatmap = ax.imshow(pd.pivot_table(series, index=series.index.time, columns=series.index.date,
                                       values="UniversalThermalClimateIndex").values[::-1],
                        extent=[dates.date2num(series.index.min()), dates.date2num(series.index.max()), 726449, 726450],
                        aspect='auto', cmap=cmap, interpolation='none', vmin=-40, vmax=46)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    ax.yaxis_date()
    ax.yaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax.invert_yaxis()
    ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(ax.get_xticklabels(), ha='left', color='#555555')
    plt.setp(ax.get_yticklabels(), color='#555555')
    [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]
    ax.grid(b=True, which='major', color='white', linestyle='--', alpha=0.9)
    cb = fig.colorbar(heatmap, cmap=cmap, norm=norm, boundaries=bounds, orientation='horizontal', drawedges=False,
                      fraction=0.05, aspect=100, pad=0.1, extend='both', ticks=[-40, -27, -13, 0, 9, 26, 32, 38, 46])
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='#555555')
    cb.outline.set_visible(False)
    ax.set_title("UTCI approximation (°C)\n{0:} - {1:} - {2:}".format(df.City[0], df.Country[0], df.StationID[0]),
                 color='#555555', y=1.03)
    ax.text(0, y_move, 'Extreme\ncold stress', ha='center', va='center', transform=ax.transAxes, color='#555555',
            fontsize='small')
    ax.text(np.interp(-27 + (-40 - -27) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\ncold stress',
            ha='center', va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(-13 + (-27 - -13) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(0 + (-13 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(0 + (9 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Slight\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(9 + (26 - 9) / 2, [-44.319, 50.319], [0, 1]), y_move, 'No thermal stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(26 + (32 - 26) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(32 + (38 - 32) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(np.interp(38 + (46 - 38) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color='#555555', fontsize='small')
    ax.text(1, y_move, 'Extreme\nheat stress', ha='center', va='center', transform=ax.transAxes, color='#555555',
            fontsize='small')
    plt.tight_layout()

    # Save figure
    if savepath != None:
        print("Saving to {}".format(savepath))
        fig.savefig(savepath, bbox_inches="tight", dpi=300, transparent=False)
    if close:
        plt.close()


def windroset(df, variable, seasonal_period, day_period, nsector=16, cmap=None, savepath=None, cls=True):
    # Descibe a set of masks to remove unwanted hours of the year
    toy_masks = {
        "Daily": ((df.index.hour >= 0) & (df.index.hour <= 24)),
        "Morning": ((df.index.hour >= 5) & (df.index.hour <= 10)),
        "Midday": ((df.index.hour >= 11) & (df.index.hour <= 13)),
        "Afternoon": ((df.index.hour >= 14) & (df.index.hour <= 18)),
        "Evening": ((df.index.hour >= 19) & (df.index.hour <= 22)),
        "Night": ((df.index.hour >= 23) | (df.index.hour <= 4)),

        "Annual": ((df.index.month >= 1) & (df.index.month <= 12)),
        "Spring": ((df.index.month >= 3) & (df.index.month <= 5)),
        "Summer": ((df.index.month >= 6) & (df.index.month <= 8)),
        "Autumn": ((df.index.month >= 9) & (df.index.month <= 11)),
        "Winter": ((df.index.month <= 2) | (df.index.month >= 12))
    }
    speed_mask = (df.wind_speed != 0)
    direction_mask = (df.wind_direction != 0)
    mask = np.array([toy_masks[day_period], toy_masks[seasonal_period], speed_mask, direction_mask]).all(axis=0)

    fig = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    if variable == "UniversalThermalClimateIndex":
        unit = "C"
        ax.bar(df.wind_direction[mask], df[variable][mask], normed=True,
               bins=[-100, -40, -27, -13, 0, 9, 26, 32, 38, 46], opening=1, edgecolor='Black',
               lw=0.25, nsector=nsector,
               colors=["#053061", "#1A5899", "#347FB9", "#82BBD9", "#BFDCEB", "#FFFFFF", "#F7C1AA", "#E3806B",
                       "#C84648", "#B2182B"])
    elif variable == "dry_bulb_temperature":
        unit = "C"
        ax.bar(df.wind_direction[mask], df[variable][mask], normed=True,
               bins=[-15, -5, 5, 15, 25, 35], opening=1, edgecolor='White',
               lw=0.25, nsector=nsector,
               cmap=plt.cm.Reds if cmap == None else cmap)
    elif variable == "relative_humidity":
        unit = "%"
        ax.bar(df.wind_direction[mask], df[variable][mask], normed=True,
               bins=[0, 20, 40, 60, 85], opening=1, edgecolor='White',
               lw=0.25, nsector=nsector,
               cmap=plt.cm.Blues if cmap == None else cmap)
    elif variable == "wind_speed":
        unit = "m/s"
        ax.bar(df.wind_direction[mask], df[variable][mask], normed=True,
               bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], opening=1, edgecolor='White',
               lw=0.25, nsector=nsector,
               cmap=plt.cm.Purples if cmap == None else cmap)

    lgd = ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', frameon=False, title=unit)
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color='#555555') for text in lgd.get_texts()]
    plt.setp(lgd.get_title(), color='#555555')

    if variable == "relative_humidity":
        lgd.get_texts()[0].set_text("Very dry (<20)")
        lgd.get_texts()[1].set_text("Dry (20-40)")
        lgd.get_texts()[2].set_text("Average (40-60)")
        lgd.get_texts()[3].set_text("Humid (60-85)")
        lgd.get_texts()[4].set_text("Very humid (>85)")
    elif variable == "UniversalThermalClimateIndex":
        lgd.get_texts()[0].set_text("Extreme cold stress")
        lgd.get_texts()[1].set_text("Very strong cold stress")
        lgd.get_texts()[2].set_text("Strong cold stress")
        lgd.get_texts()[3].set_text("Moderate cold stress")
        lgd.get_texts()[4].set_text("Slight cold stress")
        lgd.get_texts()[5].set_text("No thermal stress")
        lgd.get_texts()[6].set_text("Moderate heat stress")
        lgd.get_texts()[7].set_text("Strong heat stress")
        lgd.get_texts()[8].set_text("Very strong heat stress")
        lgd.get_texts()[9].set_text("Extreme heat stress")
    else:
        for i, leg in enumerate(lgd.get_texts()):
            b = leg.get_text().replace('[', '').replace(')', '').split(' : ')
            lgd.get_texts()[i].set_text(b[0] + ' to ' + b[1])

    ax.grid(linestyle=':', color='#555555', alpha=0.5)
    ax.spines['polar'].set_visible(False)
    plt.setp(ax.get_xticklabels(), color='#555555')
    plt.setp(ax.get_yticklabels(), color='#555555')
    ax.set_title("{2:} - {3:} - {4:}\n{0:} - {1:} - {5:}".format(df.City[0], df.Country[0], seasonal_period, day_period,
                                                                 variable, df.StationID[0]), y=1.06, color='#555555')

    plt.tight_layout()

    if savepath:
        print("Saving to {}".format(savepath))
        plt.savefig(savepath, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300, transparent=False)
    if cls:
        plt.close()