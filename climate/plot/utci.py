import matplotlib.pyplot as plt
from matplotlib import dates, cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import numpy as np

def utci_frequency_generic(series, hours=['00:00', '23:59'], tone_color="k", title=None):

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
        ax.text(i, top + 0.00025, j, ha='center', va='bottom', color=tone_color)
        ax.text(i, bottom - 0.001, k, ha='center', va='top', color=tone_color, fontsize='small')

    # Tidy plot
    plt.tight_layout()

    # Return nil data to prevent auto-plotting in Jupyter notebooks
    plt.close()

    return fig


def utci_frequency(self, variable, hours=['00:00', '23:59'], tone_color="k", save=False):

    fig = utci_frequency_generic(getattr(self, variable), hours=hours, tone_color=tone_color, title='UTCI approximation (°C) annual frequency between {2:} and {3:}\n{0:} - {1:} - {4:}'.format(self.city,
                                                                                                          self.country,
                                                                                                          hours[0],
                                                                                                          hours[1],
                                                                                                          self.station_id))

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "utci_frequency_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("UTCI frequency saved to {}".format(save_path))

    plt.close()

    return fig


def utci_heatmap_generic(series, y_move=-0.22, tone_color="k", title=None):
    sname = series.name
    series = series.to_frame()
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


def utci_heatmap(self, variable, tone_color="k", save=False):
    fig = utci_heatmap_generic(getattr(self, variable), y_move=-0.22, tone_color=tone_color,
                               title='UTCI approximation (°C)\n{0:} - {1:} - {2:}'.format(
                                     self.city,
                                     self.country,
                                     self.station_id))

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "utci_heatmap_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("UTCI heatmap saved to {}".format(save_path))

    plt.close()

    return fig
