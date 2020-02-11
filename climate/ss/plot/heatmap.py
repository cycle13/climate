from climate import renamer

import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def heatmap_generic(series, cmap='Greys', tone_color="k", title=None, vrange=None):

    # Instantiate figure
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))

    # Load data and remove timezone from index
    series = series.to_frame()

    # Reshape data into time/day matrix
    ll = series.pivot_table(columns=series.index.date, index=series.index.time).values[::-1]

    # Plot data
    heatmap = ax.imshow(
        ll,
        extent=[mdates.date2num(series.index.min()), mdates.date2num(series.index.max()), 726449, 726450],
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

    # Return nil data to prevent auto-plotting in Jupyter notebooks
    plt.close()

    return fig

def heatmap(self, variable, cmap='Greys', tone_color="k", save=False):

    fig = heatmap_generic(getattr(self, variable), cmap=cmap, tone_color=tone_color, title="{0:}\n{1:} - {2:} - {3:}".format(renamer(getattr(self, variable).name), self.city, self.country, self.station_id))

    save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "heatmap_{}.png".format(variable)

    if save:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("Heatmap saved to {}".format(save_path))

    plt.close()

    return fig
