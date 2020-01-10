import matplotlib.dates as mdates
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

months = mdates.MonthLocator()
months_fmt = mdates.DateFormatter('%b')
weeks = mdates.WeekdayLocator()


def compare_series_generic(serieses, names=None, colors=None, bins=None, tone_color="k", savepath=None, title=None):
    if names is None:
        names = ["{}_{}".format(n, i.name) for n, i in enumerate(serieses)]

    temp = pd.concat(serieses, axis=1, keys=names)

    if bins is None:
        bins = np.linspace(np.floor(temp.min().min()), np.ceil(temp.max().max()), 25)

    fig, ax = plt.subplots(2, 1, figsize=(15, 7))

    # Frequency plot
    if colors is None:
        ax[1].hist(temp.values, bins, label=temp.columns, rwidth=0.9, lw=0, density=True, zorder=4)
    else:
        ax[1].hist(temp.values, bins, label=temp.columns, color=colors, rwidth=0.9, lw=0, density=True, zorder=4)
    ax[1].set_xlabel("Value bins", color=tone_color)
    ax[1].set_ylabel("Frequency", color=tone_color)
    ax[1].set_yticklabels(['{:,.0%}'.format(x) for x in ax[1].get_yticks()], color=tone_color)
    ax[1].set_xticks(bins, minor=False)
    plt.setp(ax[1].get_xticklabels(), color=tone_color)
    ax[1].grid(axis="both", which="major", color=tone_color, ls="--", alpha=0.2, zorder=2)
    ax[1].set_xlim(np.floor(temp.min().min()), np.ceil(temp.max().max()))

    # Series plot
    if colors is None:
        ax[0].plot(temp, lw=1, zorder=4)
    else:
        for n, i in enumerate(serieses):
            ax[0].plot(temp.index, i, lw=1, zorder=4, color=colors[n])
    plt.setp(ax[0].get_xticklabels(), ha='left', color=tone_color)
    plt.setp(ax[0].get_yticklabels(), color=tone_color)
    ax[0].set_ylabel("Value", color=tone_color)
    ax[0].set_xlim(temp.index.min(), temp.index.max())
    [ax[0].spines[j].set_visible(False) for j in ["right", "top"]]
    ax[0].xaxis.set_major_locator(months)
    ax[0].xaxis.set_major_formatter(months_fmt)
    ax[0].xaxis.set_minor_locator(weeks)
    ax[0].grid(axis="both", which="major", color=tone_color, ls="--", alpha=0.2, zorder=2)

    [i.tick_params(length=0) for i in ax]
    [[i.spines[j].set_visible(False) for j in ["right", "top"]] for i in ax]
    [[i.spines[j].set_color(tone_color) for j in ["left", "bottom"]] for i in ax]

    lgd = ax[1].legend(
        bbox_to_anchor=(0, 1),
        loc=2,
        ncol=1,
        borderaxespad=0,
        frameon=False,
    )
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color=tone_color) for text in lgd.get_texts()]

    # Add title if provided
    if title is not None:
        plt.suptitle(title, color=tone_color, y=1)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=300)

    plt.close()

    return fig