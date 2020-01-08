import matplotlib.pyplot as plt
from matplotlib import dates, cm
from matplotlib.ticker import MaxNLocator
from climate.common.constants import utci_cmap
from matplotlib.colors import ListedColormap, BoundaryNorm

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import pandas as pd
import numpy as np

def utci_frequency_generic(series, hours=['00:00', '23:59'], tone_color="k", title=None, cmap=None):

    if cmap is None:
        cmap = utci_cmap
    else:
        try:
            cmap = cm.get_cmap(cmap)
        except Exception as e:
            print(e)
            cmap = cmap

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    data = series.between_time(hours[0], hours[1], include_start=True, include_end=False)
    data.plot.hist(bins=np.arange(-40, 58, 1), ax=ax, zorder=4, color="#555555", alpha=0.8, density=True)
    ax.set_yticklabels(['{:0.1f}%'.format(x * 100) for x in ax.get_yticks()])
    ax.set_xlabel('UTCI (C)', color=tone_color, labelpad=20)
    ax.tick_params(axis='both', colors=tone_color)
    ax.set_ylabel('Frequency', color=tone_color)
    ax.set_xticks([-40, -27, -13, 0, 9, 26, 32, 38, 46])
    ax.tick_params(axis='both', which='major')
    ax.grid(b=True, which='major', color='white', linestyle='--', alpha=0.9, zorder=10)
    [ax.spines[spine].set_visible(False) for spine in ['top', 'right']]
    [ax.spines[j].set_color(tone_color) for j in ['bottom', 'left']]

    # Add title if provided
    if title is not None:
        plt.title(title, color=tone_color, y=1.05)

    ax.axvspan(-50, -40, color=cmap(-0.1), zorder=1)
    ax.axvspan(-40, -27, color=cmap(0), zorder=1)
    ax.axvspan(-27, -13, color=cmap(0.2), zorder=1)
    ax.axvspan(-13, 0, color=cmap(0.4), zorder=1)
    ax.axvspan(0, 9, color=cmap(0.5), zorder=1)
    ax.axvspan(9, 26, color=cmap(0.7), zorder=1)
    ax.axvspan(26, 32, color=cmap(0.8), zorder=1)
    ax.axvspan(32, 38, color=cmap(0.9), zorder=1)
    ax.axvspan(38, 46, color=cmap(0.95), zorder=1)
    ax.axvspan(46, 60, color=cmap(1.2), zorder=1)
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
        ax.text(i, top + 0.00025, j, ha='center', va='bottom', color=tone_color)
        ax.text(i, bottom - 0.001, k, ha='center', va='top', color=tone_color, fontsize='small')

    # Tidy plot
    plt.tight_layout()

    # Return nil data to prevent auto-plotting in Jupyter notebooks
    plt.close()

    return fig


def utci_frequency(self, variable, hours=['00:00', '23:59'], tone_color="k", save=False, cmap=None):

    fig = utci_frequency_generic(getattr(self, variable), hours=hours, tone_color=tone_color, title='UTCI approximation (°C) annual frequency between {2:} and {3:}\n{0:} - {1:} - {4:}'.format(self.city,
                                                                                                          self.country,
                                                                                                          hours[0],
                                                                                                          hours[1],
                                                                                                          self.station_id), cmap=cmap)

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "utci_frequency_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("UTCI frequency saved to {}".format(save_path))

    plt.close()

    return fig


def utci_heatmap_generic(series, y_move=-0.22, tone_color="k", title=None, cmap=None):
    sname = series.name
    series = series.to_frame()

    if cmap is None:
        cmap = utci_cmap
    else:
        try:
            cmap = cm.get_cmap(cmap)
        except Exception as e:
            cmap = cmap

    bounds = np.arange(-41, 48, 1)
    norm = BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    heatmap = ax.imshow(pd.pivot_table(series, index=series.index.time, columns=series.index.date,
                                       values=sname).values[::-1],
                        extent=[dates.date2num(series.index.min()), dates.date2num(series.index.max()), 726449, 726450],
                        aspect='auto', cmap=cmap, interpolation='none', vmin=-40, vmax=46)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    ax.yaxis_date()
    ax.yaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax.invert_yaxis()
    ax.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(ax.get_xticklabels(), ha='left', color=tone_color)
    plt.setp(ax.get_yticklabels(), color=tone_color)
    [ax.spines[spine].set_visible(False) for spine in ['top', 'bottom', 'left', 'right']]
    ax.grid(b=True, which='major', color='white', linestyle='--', alpha=0.9)
    cb = fig.colorbar(heatmap, cmap=cmap, norm=norm, boundaries=bounds, orientation='horizontal', drawedges=False,
                      fraction=0.05, aspect=100, pad=0.1, extend='both', ticks=[-40, -27, -13, 0, 9, 26, 32, 38, 46])
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=tone_color)
    cb.outline.set_visible(False)

    # Add title if provided
    if title is not None:
        plt.title(title, color=tone_color, y=1.03)

    ax.text(0, y_move, 'Extreme\ncold stress', ha='center', va='center', transform=ax.transAxes, color=tone_color,
            fontsize='small')
    ax.text(np.interp(-27 + (-40 - -27) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\ncold stress',
            ha='center', va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(-13 + (-27 - -13) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(0 + (-13 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(0 + (9 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Slight\ncold stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(9 + (26 - 9) / 2, [-44.319, 50.319], [0, 1]), y_move, 'No thermal stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(26 + (32 - 26) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(32 + (38 - 32) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(np.interp(38 + (46 - 38) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\nheat stress', ha='center',
            va='center', transform=ax.transAxes, color=tone_color, fontsize='small')
    ax.text(1, y_move, 'Extreme\nheat stress', ha='center', va='center', transform=ax.transAxes, color=tone_color,
            fontsize='small')
    plt.tight_layout()

    return fig


def utci_heatmap(self, variable, tone_color="k", save=False, cmap=None):

    fig = utci_heatmap_generic(getattr(self, variable), y_move=-0.22, tone_color=tone_color,
                               title='UTCI approximation (°C)\n{0:} - {1:} - {2:}'.format(
                                     self.city,
                                     self.country,
                                     self.station_id), cmap=cmap)

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "utci_heatmap_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("UTCI heatmap saved to {}".format(save_path))

    plt.close()

    return fig


def utci_heatmap_detailed_generic(series, tone_color="k", title=None, cmap=None):
    sname = series.name
    series = series.to_frame()

    if cmap is None:
        cmap = utci_cmap
    else:
        try:
            cmap = cm.get_cmap(cmap)
        except Exception as e:
            print(e)
            cmap = cmap

    fig = plt.figure(figsize=(15, 6), constrained_layout=True)
    spec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[2, 1], hspace=0.1)

    hmap = fig.add_subplot(spec[0, 0])
    hbar = fig.add_subplot(spec[1, 0])

    bounds = np.arange(-41, 48, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    # Plot heatmap
    heatmap = hmap.imshow(pd.pivot_table(series, index=series.index.time, columns=series.index.date,
                                         values=sname).values[::-1],
                          extent=[dates.date2num(series.index.min()), dates.date2num(series.index.max()), 726449,
                                  726450],
                          aspect='auto', cmap=cmap, interpolation='none', vmin=-40, vmax=46)
    hmap.xaxis_date()
    hmap.xaxis.set_major_formatter(dates.DateFormatter('%b'))
    hmap.yaxis_date()
    hmap.yaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    hmap.invert_yaxis()
    hmap.tick_params(labelleft=True, labelright=True, labelbottom=True)
    plt.setp(hmap.get_xticklabels(), ha='left', color=tone_color)
    plt.setp(hmap.get_yticklabels(), color=tone_color)
    for spine in ['top', 'bottom', 'left', 'right']:
        hmap.spines[spine].set_visible(False)
        hmap.spines[spine].set_color(tone_color)
    hmap.grid(b=True, which='major', color=tone_color, linestyle=':', alpha=0.5)

    # Add colorbar legend and text descriptors for comfort bands
    cb = fig.colorbar(heatmap, cmap=cmap, norm=norm, boundaries=bounds,
                      orientation='horizontal', drawedges=False, fraction=0.01, aspect=50,
                      pad=-0.0, extend='both', ticks=[-40, -27, -13, 0, 9, 26, 32, 38, 46])
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=tone_color)
    cb.outline.set_visible(False)
    # cb.outline.set_color("#555555")
    y_move = -0.4
    hbar.text(0, y_move, 'Extreme\ncold stress', ha='center', va='center', transform=hbar.transAxes, color=tone_color,
              fontsize='small')
    hbar.text(np.interp(-27 + (-40 - -27) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\ncold stress',
              ha='center', va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(-13 + (-27 - -13) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\ncold stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(0 + (-13 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\ncold stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(0 + (9 - 0) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Slight\ncold stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(9 + (26 - 9) / 2, [-44.319, 50.319], [0, 1]), y_move, 'No thermal stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(26 + (32 - 26) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Moderate\nheat stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(32 + (38 - 32) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Strong\nheat stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(np.interp(38 + (46 - 38) / 2, [-44.319, 50.319], [0, 1]), y_move, 'Very strong\nheat stress', ha='center',
              va='center', transform=hbar.transAxes, color=tone_color, fontsize='small')
    hbar.text(1, y_move, 'Extreme\nheat stress', ha='center', va='center', transform=hbar.transAxes, color='#555555',
              fontsize='small')

    # Add stacked plot
    bins = [-100, -40, -27, -13, 0, 9, 26, 32, 38, 46, 100]
    tags = ["Extreme cold stress", "Very strong cold stress", "Strong cold stress", "Moderate cold stress",
            "Slight cold stress", "No thermal stress", "Moderate heat stress", "Strong heat stress",
            "Very strong heat stress", "Extreme heat stress"]
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


    clrs = [
        cmap(-0.1),
        cmap(0),
        cmap(0.2),
        cmap(0.4),
        cmap(0.5),
        cmap(0.7),
        cmap(0.8),
        cmap(0.9),
        cmap(0.95),
        cmap(1.2),
    ]

    adf = pd.DataFrame()
    for mnth_n, mnth in enumerate(months):
        # Filter the series to return only the month
        a = series[series.index.month == mnth_n + 1].dropna().values
        a = pd.Series(index=tags, name=mnth,
                      data=[((a > i) & (a <= j)).sum() / len(a) for n, (i, j) in enumerate(zip(bins[:-1], bins[1:]))])
        adf = pd.concat([adf, a], axis=1)
    adf = adf.T[tags]
    adf.plot(kind="bar", ax=hbar, stacked=True, colors=clrs, width=1, legend=False)
    hbar.set_xlim(-0.5, 11.5)
    plt.setp(hbar.get_xticklabels(), ha='center', rotation=0, color=tone_color)
    plt.setp(hbar.get_xticklabels(), ha='left', color=tone_color)
    plt.setp(hbar.get_yticklabels(), color=tone_color)
    for spine in ['top', 'right']:
        hbar.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        hbar.spines[spine].set_color(tone_color)
    hbar.grid(b=True, which='major', color=tone_color, linestyle=':', alpha=0.5)
    hbar.set_yticklabels(['{:,.0%}'.format(x) for x in hbar.get_yticks()])

    # Add header percentages for bar plot
    cold_percentages = adf.iloc[:, :5].sum(axis=1).values
    comfortable_percentages = adf.iloc[:, 5]
    hot_percentages = adf.iloc[:, 6:].sum(axis=1).values
    for n, (i, j, k) in enumerate(zip(*[cold_percentages, comfortable_percentages, hot_percentages])):
        hbar.text(n, 1.02, "{0:0.1f}%\n\n".format(i * 100), va="bottom", ha="center", color=cmap(-0.1), fontsize="small")
        hbar.text(n, 1.02, "{0:0.1f}%\n".format(j * 100), va="bottom", ha="center", color="#555555", fontsize="small")
        hbar.text(n, 1.02, "{0:0.1f}%".format(k * 100), va="bottom", ha="center", color=cmap(1.2), fontsize="small")
    hbar.set_ylim(0, 1)

    # Add title if provided
    if title is not None:
        plt.suptitle(title, color=tone_color, y=1, va="bottom")

    plt.tight_layout()

    return fig


def utci_heatmap_detailed(self, variable, tone_color="k", save=False, cmap=cmap):
    fig = utci_heatmap_detailed_generic(getattr(self, variable), tone_color=tone_color,
                               title='UTCI approximation (°C)\n{0:} - {1:} - {2:}'.format(
                                     self.city,
                                     self.country,
                                     self.station_id), cmap=cmap)

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "utci_heatmap_detailed_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("UTCI heatmap (detailed) saved to {}".format(save_path))

    plt.close()

    return fig
