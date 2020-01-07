import matplotlib.pyplot as plt
from matplotlib import dates, cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
from windrose import WindroseAxes
import pandas as pd
import numpy as np

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


def windrose(epw_object, seasonal_period, day_period, nsector=16, cmap=None, savepath=None, cls=True):
    # Descibe a set of masks to remove unwanted hours of the year
    toy_masks = {
        "Daily": ((epw_object.df.index.hour >= 0) & (epw_object.df.index.hour <= 24)),
        "Morning": ((epw_object.df.index.hour >= 5) & (epw_object.df.index.hour <= 10)),
        "Midday": ((epw_object.df.index.hour >= 11) & (epw_object.df.index.hour <= 13)),
        "Afternoon": ((epw_object.df.index.hour >= 14) & (epw_object.df.index.hour <= 18)),
        "Evening": ((epw_object.df.index.hour >= 19) & (epw_object.df.index.hour <= 22)),
        "Night": ((epw_object.df.index.hour >= 23) | (epw_object.df.index.hour <= 4)),

        "Annual": ((epw_object.df.index.month >= 1) & (epw_object.df.index.month <= 12)),
        "Spring": ((epw_object.df.index.month >= 3) & (epw_object.df.index.month <= 5)),
        "Summer": ((epw_object.df.index.month >= 6) & (epw_object.df.index.month <= 8)),
        "Autumn": ((epw_object.df.index.month >= 9) & (epw_object.df.index.month <= 11)),
        "Winter": ((epw_object.df.index.month <= 2) | (epw_object.df.index.month >= 12))
    }
    speed_mask = (epw_object.wind_speed != 0)
    direction_mask = (epw_object.wind_direction != 0)
    mask = np.array([toy_masks[day_period], toy_masks[seasonal_period], speed_mask, direction_mask]).all(axis=0)

    fig = plt.figure(figsize=(6, 6))
    ax = WindroseAxes.from_ax()
    unit = "m/s"
    ax.bar(epw_object.wind_direction[mask], epw_object.wind_speed[mask], normed=True,
           bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], opening=1, edgecolor='White',
           lw=0.25, nsector=nsector,
           cmap=plt.cm.Purples if cmap is None else cmap)

    lgd = ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', frameon=False, title=unit)
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color='#555555') for text in lgd.get_texts()]
    plt.setp(lgd.get_title(), color='#555555')

    for i, leg in enumerate(lgd.get_texts()):
        b = leg.get_text().replace('[', '').replace(')', '').split(' : ')
        lgd.get_texts()[i].set_text(b[0] + ' to ' + b[1])

    ax.grid(linestyle=':', color='#555555', alpha=0.5)
    ax.spines['polar'].set_visible(False)
    plt.setp(ax.get_xticklabels(), color='#555555')
    plt.setp(ax.get_yticklabels(), color='#555555')
    ax.set_title("{2:} - {3:} - {4:}\n{0:} - {1:} - {5:}".format(epw_object.city, epw_object.country, seasonal_period, day_period, "Wind speed", epw_object.station_id), y=1.06, color='#555555')

    plt.tight_layout()

    if savepath:
        print("Saving to {}".format(savepath))
        plt.savefig(savepath, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300, transparent=True)
    if cls:
        plt.close()