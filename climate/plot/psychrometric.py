import matplotlib.pyplot as plt
import numpy as np
from climate.compute.psychrometrics import humidity_ratio_from_vapor_pressure, humidity_ratio_from_relative_humidity, humidity_ratio_from_enthalpy, moist_air_enthalpy, moist_air_volume, dew_point_from_vapor_pressue, dry_bulb_temperature_from_enthalpy, dry_bulb_temperature_from_specific_volume, wet_bulb_temperature_from_humidity_ratio, saturation_vapor_pressure, degree_of_saturation, vapor_pressure_from_humidity_ratio

# TODO: Use methods from psychrochart for plotting of enthalpy, specific volume and wet-bulb curves


def psychrometric(self, bins=50, cmap="Greys", tone_color="k", save=False):
    # Check that humidity ratio is available - if not then throw error
    if self.humidity_ratio is None:
        raise Exception('Psychrometric calculations have not been run yet. Try running self.annual_psychrometrics()')

    # Construct the save_path and create directory if it doesn't exist
    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "psychrometrics.png"

    ## Internal fixed variables
    pressure = 101325

    dbt_limits = [-20, 50]
    dbt_step = 5
    dbt_detailed_step = 1
    dbt_range = np.arange(dbt_limits[0], dbt_limits[1] + dbt_step, dbt_step)
    dbt_range_detailed = np.arange(dbt_limits[0], dbt_limits[1] + dbt_detailed_step, dbt_detailed_step)

    hr_range = np.arange(0, 0.035, 0.005)
    rh_range = np.arange(0, 110, 10)
    enthalpy_range = np.arange(5, 105, 5)

    # Instantiate plot
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))

    # Canvas formatting
    [ax.spines[spine].set_visible(False) for spine in ['top', 'left']]
    [ax.spines[spine].set_color(tone_color) for spine in ['bottom', 'right']]

    # Plot constant lines for Dry-Bulb Temperature
    saturation_humidity_ratio = humidity_ratio_from_vapor_pressure(saturation_vapor_pressure(dbt_range), pressure)
    min_humidity_ratio = np.repeat(hr_range[0], len(saturation_humidity_ratio))
    ax.plot([dbt_range, dbt_range], [min_humidity_ratio, saturation_humidity_ratio], color=tone_color, alpha=0.5, lw=0.5, zorder=8)

    # Plot constant lines for Humidity Ratio
    min_x = dew_point_from_vapor_pressue(dbt_range[-1], vapor_pressure_from_humidity_ratio(hr_range, pressure))
    max_x = np.repeat(dbt_range[-1], len(min_x))
    ax.plot([min_x, max_x], [hr_range, hr_range], color=tone_color, alpha=0.5, lw=0.5, zorder=8)

    # Plot constant lines for Relative Humidity
    rh_ys = np.array([humidity_ratio_from_relative_humidity(dbt_range_detailed, i / 100, pressure) for i in rh_range])
    for i in range(len(rh_ys)):
        ax.plot(np.tile(dbt_range_detailed, (len(rh_ys), 1))[i], rh_ys[i], color=tone_color, alpha=0.5, lw=0.5, zorder=8)

    # TODO - Plot constant lines for Enthalpy


    # TODO - Plot constant lines for Wet-Bulb Temperature

    # TODO - Plot constant lines for Specific Volume


    # Plot values
    counts, xedges, yedges, im = ax.hist2d(
        self.dry_bulb_temperature, self.humidity_ratio, bins=bins, cmin=1, alpha=0.9, normed=False, cmap=cmap, lw=0, zorder=1)

    # Fill top part of chart to tidy edge of values
    ax.fill_between(dbt_range_detailed, rh_ys[-1], 100, interpolate=True, color='w', lw=0, edgecolor=None, zorder=4)

    # Colorbar
    cb = plt.colorbar(im, ax=ax, shrink=1, pad=0.071)
    cb.ax.set_title('Hours', color=tone_color)
    cb.outline.set_visible(False)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tone_color)

    # X axis formatting
    ax.set_xlim(dbt_range[0], dbt_range[1])
    ax.set_xticks(dbt_range)
    ax.set_xlabel("Dry-bulb temperature (°C)", color=tone_color, fontsize="x-large")
    plt.setp(plt.getp(ax.axes, 'xticklabels'), color=tone_color)

    # Y axis formatting
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_ylim(hr_range[0], hr_range[1])
    ax.set_yticks([i for i in hr_range])
    ax.set_ylabel("Humidity ratio ($kg_{water}/kg_{air}$)", color=tone_color, fontsize="x-large")
    plt.setp(plt.getp(ax.axes, 'yticklabels'), color=tone_color)

    # Generate summary metrics
    summary_text_size = "medium"
    summary_date_format = "%b %d %H:%M"
    summary_spacing = {"border": 0.0, "h": 0.22, "v": 0.24, "top": 0.02}

    # Peak cooling
    max_dbt_idx = self.dry_bulb_temperature.idxmax()
    max_dbt_table = "Peak cooling {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(max_dbt_idx.strftime(summary_date_format), self.wind_speed[max_dbt_idx], self.wind_direction[max_dbt_idx], self.dry_bulb_temperature[max_dbt_idx], self.wet_bulb_temperature[max_dbt_idx], self.relative_humidity[max_dbt_idx], self.dew_point_temperature[max_dbt_idx], self.enthalpy[max_dbt_idx] / 1000, self.humidity_ratio[max_dbt_idx])
    ax.text(0 + summary_spacing["border"], 1 - summary_spacing["top"] - summary_spacing["border"], max_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=summary_text_size, color=tone_color, **{'fontname': 'monospace'})

    # Peak heating
    min_dbt_idx = self.dry_bulb_temperature.idxmin()
    min_dbt_table = "Peak heating {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(min_dbt_idx.strftime(summary_date_format), self.wind_speed[min_dbt_idx], self.wind_direction[min_dbt_idx], self.dry_bulb_temperature[min_dbt_idx], self.wet_bulb_temperature[min_dbt_idx], self.relative_humidity[min_dbt_idx], self.dew_point_temperature[min_dbt_idx], self.enthalpy[min_dbt_idx] / 1000, self.humidity_ratio[min_dbt_idx])
    ax.text(0 + summary_spacing["border"], 1 - summary_spacing["top"] - summary_spacing["border"] - summary_spacing["v"], min_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=summary_text_size, color=tone_color, **{'fontname': 'monospace'})

    # Peak humidity-ratio
    max_hr_idx = self.humidity_ratio.idxmax()
    max_hr_table = "Peak humidity ratio {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(max_hr_idx.strftime(summary_date_format), self.wind_speed[max_hr_idx], self.wind_direction[max_hr_idx], self.dry_bulb_temperature[max_hr_idx], self.wet_bulb_temperature[max_hr_idx], self.relative_humidity[max_hr_idx], self.dew_point_temperature[max_hr_idx], self.enthalpy[max_hr_idx] / 1000, self.humidity_ratio[max_hr_idx])
    ax.text(summary_spacing["h"], 1 - summary_spacing["top"] - summary_spacing["border"], max_hr_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=summary_text_size, color=tone_color, **{'fontname': 'monospace'})

    # Generate max enthalpy summary
    max_enthalpy_idx = self.enthalpy.idxmax()
    max_enthalpy_table = "Peak enthalpy ratio {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(max_enthalpy_idx.strftime(summary_date_format), self.wind_speed[max_enthalpy_idx], self.wind_direction[max_enthalpy_idx], self.dry_bulb_temperature[max_enthalpy_idx], self.wet_bulb_temperature[max_enthalpy_idx], self.relative_humidity[max_enthalpy_idx], self.dew_point_temperature[max_enthalpy_idx], self.enthalpy[max_enthalpy_idx] / 1000, self.humidity_ratio[max_enthalpy_idx])
    ax.text(summary_spacing["h"], 1 - summary_spacing["top"] - summary_spacing["border"] - summary_spacing["v"], max_enthalpy_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=summary_text_size, color=tone_color, **{'fontname': 'monospace'})

    # Add title
    ti = ax.set_title("{} - {} - {}".format(self.city, self.country, self.station_id), color=tone_color, loc="left", fontsize="xx-large")

    # Add acronym keys
    keys = "WS: Wind speed | WD: Wind direction | DBT: Dry-bulb temperature | WBT: Wet-bulb temperature\nRH: Relative humidity | DPT: Dew-point temperature | h: Enthalpy | HR: Humidity ratio"
    te = ax.text(0.5, -0.1, keys, transform=ax.transAxes, ha="center", va="top", zorder=8, fontsize="medium", color=tone_color, **{'fontname': 'monospace'})

    plt.tight_layout()

    # Save figure
    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False, bbox_extra_artists=[ti, te, ])
        print("Psychrometric plot saved to {}".format(save_path))

    plt.close()

    return fig

#
# def psychrometric_(self, bins=50, cmap="Greys", tone_color="k", save=False):
#
#     # Check that humidity ratio is available - if not then throw error
#     if self.humidity_ratio is None:
#         raise Exception('Psychrometric calculations have not been run yet. Try running self.annual_psychrometrics()')
#
#     # Construct the save_path and create directory if it doesn't exist
#     save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "psychrometrics.png"
#
#     ## Internal fixed variables
#     pressure = 101325
#
#     dbt_limits = [-20, 50]
#     dbt_step = 5
#     dbt_detailed_step = 1
#     dbt_range = np.arange(dbt_limits[0], dbt_limits[1] + dbt_step, dbt_step)
#     dbt_range_detailed = np.arange(dbt_limits[0], dbt_limits[1] + dbt_detailed_step, dbt_detailed_step)
#
#     hr_range = np.arange(0, 0.035, 0.005)
#     rh_range = np.arange(0, 110, 10)
#     #
#     # dry_bulb_temperature_plot_range = range(-20, 51, 1)
#     # humidity_ratio_plot_range = [i / 10000 for i in range(0, 301, 1)]
#     # enthalpy_plot_range = range(-10, 120, 10)
#     # relative_humidity_plot_range = [i / 100 for i in range(0, 101, 10)]
#     #
#     # # figure instantiation
#     # fig, ax = plt.subplots(1, 1, figsize=(15, 8))
#     #
#     # # plot values from weather file
#     # counts, xedges, yedges, im = ax.hist2d(self.dry_bulb_temperature, self.humidity_ratio, bins=bins, cmin=1, alpha=0.9, normed=False, cmap=cmap, lw=0,
#     #                                        zorder=0)
#     #
#     # # y-axis formatting
#     # ax.yaxis.tick_right()
#     # ax.yaxis.set_label_position("right")
#     # ax.set_ylim(0, 0.03)
#     # ax.set_yticks([i / 1000 for i in range(0, 35, 5)])
#     # ax.set_ylabel("Humidity ratio ($kg_{water}/kg_{air}$)", color=tone_color, fontsize="x-large")
#     #
#     # # x-axis formatting
#     # ax.set_xlim(-20, 50)
#     # ax.set_xticks(range(-20, 55, 5))
#     # ax.set_xlabel("Dry-bulb temperature ($°C$)", color=tone_color, fontsize="x-large")
#     #
#     # # canvas formatting
#     # ax.tick_params(axis="both", color=tone_color, grid_color=tone_color, grid_alpha=1, grid_lw=0.5)
#     # [ax.spines[spine].set_visible(False) for spine in ['top', 'left']]
#     # [ax.spines[spine].set_color(tone_color) for spine in ['bottom', 'right']]
#     # [ax.spines[spine].set_lw(1) for spine in ['top', 'left', 'bottom', 'right']]
#     #
#     # # relative humidity grid/curves
#     # n = 0
#     # for rh in relative_humidity_plot_range:
#     #     h_r = [GetHumRatioFromRelHum(i, rh, 101325) for i in dry_bulb_temperature_plot_range]
#     #     ax.plot(dry_bulb_temperature_plot_range, h_r, color=tone_color, alpha=1, lw=0.2)
#     #     # Fill the top part of the plot
#     #     if rh == 1:
#     #         ax.fill_between(dry_bulb_temperature_plot_range, h_r, 0.031, interpolate=True, color='w', lw=0,
#     #                         edgecolor=None,
#     #                         zorder=4)
#     #     # add annotation describing line
#     #     ax.text(30, GetHumRatioFromRelHum(30, rh, 101350) + 0.0, "{0:0.0f}% RH".format(rh * 100),
#     #             ha="right", va="bottom", rotation=0, zorder=9, fontsize="small", color=tone_color)  # n*55
#     #     n += 1 / len(relative_humidity_plot_range)
#
#     # TODO: Fix enthalpy grid curves
#     # # enthalpy grid/curves
#     # for enthalpy in enthalpy_plot_range:
#     #     ys = [0, 0.030]
#     #     xs = np.array([GetTDryBulbFromEnthalpyAndHumRatio(enthalpy, i) for i in ys]) / 1000
#     #     if (enthalpy <= 50) & (enthalpy != 30):
#     #         ax.text(xs[0], 0.0002, "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="k", zorder=9,
#     #                 fontsize="small")
#     #     else:
#     #         pass
#     #     # ax.text(50, ys[0], "{}kJ/kg".format(enthalpy), ha="right", va="bottom", color="#555555", zorder=9, fontsize="small")
#     #     ax.plot(xs, ys, color="k", alpha=1, lw=0.2)
#
#     # # grid formatting
#     # ax.grid(True, lw=0.2, zorder=5, color=tone_color)
#     #
#     # # Generate summary metrics
#     # text_fontsize = "medium"
#     # date_format = "%m/%d %H:%M"
#     #
#     # # Generate peak cooling summary
#     # max_dbt_idx = self.dry_bulb_temperature.idxmax()
#     # max_dbt_table = "Peak cooling {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(max_dbt_idx.strftime(date_format), self.wind_speed[max_dbt_idx], self.wind_direction[max_dbt_idx], self.dry_bulb_temperature[max_dbt_idx], self.wet_bulb_temperature[max_dbt_idx], self.relative_humidity[max_dbt_idx], self.dew_point_temperature[max_dbt_idx], self.enthalpy[max_dbt_idx] / 1000, self.humidity_ratio[max_dbt_idx])
#     # ax.text(0, 0.98, max_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize, color=tone_color, **{'fontname': 'monospace'})
#     #
#     # # Generate peak heating summary
#     # min_dbt_idx = self.dry_bulb_temperature.idxmin()
#     # min_dbt_table = "Peak heating {0:}\nWS:  {1:>6.1f} m/s\nWD:  {2:>6.1f} deg\nDBT: {3:>6.1f} °C\nWBT: {4:>6.1f} °C\nRH:  {5:>6.1f} %\nDPT: {6:>6.1f} °C\nh:   {7:>6.1f} kJ/kg\nHR:  {8:<5.4f} kg/kg".format(
#     #     min_dbt_idx.strftime(date_format),
#     #     self.wind_speed[min_dbt_idx],
#     #     self.wind_direction[min_dbt_idx],
#     #     self.dry_bulb_temperature[min_dbt_idx],
#     #     self.wet_bulb_temperature[min_dbt_idx],
#     #     self.relative_humidity[min_dbt_idx],
#     #     self.dew_point_temperature[min_dbt_idx],
#     #     self.enthalpy[min_dbt_idx] / 1000,
#     #     self.humidity_ratio[min_dbt_idx])
#     #
#     # ax.text(0, 0.72, min_dbt_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize, color=tone_color, **{'fontname': 'monospace'})
#     #
#     # # Generate max HumidityRatio summary
#     # max_hr_idx = self.humidity_ratio.idxmax()
#     # max_hr_table = "Peak humidity ratio {0:}\n" \
#     #                "WS:  {1:>6.1f} m/s\n" \
#     #                "WD:  {2:>6.1f} deg\n" \
#     #                "DBT: {3:>6.1f} °C\n" \
#     #                "WBT: {4:>6.1f} °C\n" \
#     #                "RH:  {5:>6.1f} %\n" \
#     #                "DPT: {6:>6.1f} °C\n" \
#     #                "h:   {7:>6.1f} kJ/kg\n" \
#     #                "HR:  {8:<5.4f} kg/kg".format(
#     #     max_hr_idx.strftime(date_format),
#     #     self.wind_speed[max_hr_idx],
#     #     self.wind_direction[max_hr_idx],
#     #     self.dry_bulb_temperature[max_hr_idx],
#     #     self.wet_bulb_temperature[max_hr_idx],
#     #     self.relative_humidity[max_hr_idx],
#     #     self.dew_point_temperature[max_hr_idx],
#     #     self.enthalpy[max_hr_idx] / 1000,
#     #     self.humidity_ratio[max_hr_idx])
#     #
#     # ax.text(0.20, 0.98, max_hr_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize, color=tone_color, **{'fontname': 'monospace'})
#     #
#     # # Generate max enthalpy summary
#     # max_enthalpy_idx = self.enthalpy.idxmax()
#     # max_enthalpy_table = "Peak enthalpy ratio {0:}\n" \
#     #                      "WS:  {1:>6.1f} m/s\n" \
#     #                      "WD:  {2:>6.1f} deg\n" \
#     #                      "DBT: {3:>6.1f} °C\n" \
#     #                      "WBT: {4:>6.1f} °C\n" \
#     #                      "RH:  {5:>6.1f} %\n" \
#     #                      "DPT: {6:>6.1f} °C\n" \
#     #                      "h:   {7:>6.1f} kJ/kg\n" \
#     #                      "HR:  {8:<5.4f} kg/kg".format(
#     #     max_enthalpy_idx.strftime(date_format),
#     #     self.wind_speed[max_enthalpy_idx],
#     #     self.wind_direction[max_enthalpy_idx],
#     #     self.dry_bulb_temperature[max_enthalpy_idx],
#     #     self.wet_bulb_temperature[max_enthalpy_idx],
#     #     self.relative_humidity[max_enthalpy_idx],
#     #     self.dew_point_temperature[max_enthalpy_idx],
#     #     self.enthalpy[max_enthalpy_idx] / 1000,
#     #     self.humidity_ratio[max_enthalpy_idx])
#     # ax.text(0.20, 0.72, max_enthalpy_table, transform=ax.transAxes, ha="left", va="top", zorder=8, fontsize=text_fontsize, color=tone_color, **{'fontname': 'monospace'})
#     #
#     # # Title formatting
#     # ti = ax.set_title("{} - {} - {}".format(self.city, self.country, self.station_id),
#     #                   color=tone_color, loc="left", fontsize="xx-large")
#     #
#     # # text
#     # keys = "WS: Wind speed | WD: Wind direction | DBT: Dry-bulb temperature | WBT: Wet-bulb temperature\nRH: Relative humidity | DPT: Dew-point temperature | h: Enthalpy | HR: Humidity ratio"
#     # te = ax.text(0.5, -0.1, keys, transform=ax.transAxes, ha="center", va="top", zorder=8, fontsize="medium",
#     #              color=tone_color, **{'fontname': 'monospace'})
#     #
#     # # Colorbar
#     # cb = plt.colorbar(im, ax=ax, shrink=1, pad=0.071)
#     # cb.ax.set_title('Hours', color=tone_color)
#     # cb.outline.set_visible(False)
#     # plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=tone_color)
#     #
#     # plt.tight_layout()
#     #
#     # # Save figure
#     # if save:
#     #     save_path.parent.mkdir(parents=True, exist_ok=True)
#     #     fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False, bbox_extra_artists=[ti, te, ])
#     #     print("Psychrometric plot saved to {}".format(save_path))
#
#     return fig
