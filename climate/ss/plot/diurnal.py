from climate import renamer

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

def diurnal(self, dew_point=True, tone_color="k", save=False):

    grouping = [self.index.month, self.index.hour]

    # Group dry-bulb temperatures
    dbt_grouping = self.dry_bulb_temperature.groupby(grouping)
    dbt_min = dbt_grouping.min().reset_index(drop=True)
    dbt_mean = dbt_grouping.mean().reset_index(drop=True)
    dbt_max = dbt_grouping.max().reset_index(drop=True)
    dbt_name = renamer(self.dry_bulb_temperature.name)

    # Group relative humidity / dew-point temperature
    if dew_point:
        moisture_grouping = self.dew_point_temperature.groupby(grouping)
        moisture_name = renamer(self.dew_point_temperature.name)
    else:
        moisture_grouping = self.relative_humidity.groupby(grouping)
        moisture_name = renamer(self.relative_humidity.name)

    moisture_min = moisture_grouping.min().reset_index(drop=True)
    moisture_mean = moisture_grouping.mean().reset_index(drop=True)
    moisture_max = moisture_grouping.max().reset_index(drop=True)

    # Group solar radiation
    global_radiation_mean = self.global_horizontal_radiation.groupby(grouping).mean().reset_index(drop=True)
    glob_rad_name = renamer(self.global_horizontal_radiation.name)
    diffuse_radiation_mean = self.diffuse_horizontal_radiation.groupby(grouping).mean().reset_index(drop=True)
    dif_rad_name = renamer(self.diffuse_horizontal_radiation.name)
    direct_radiation_mean = self.direct_normal_radiation.groupby(grouping).mean().reset_index(drop=True)
    dir_rad_name = renamer(self.direct_normal_radiation.name)

    # Instantiate plot
    fig, ax = plt.subplots(3, 1, figsize=(15, 8))

    # Plot DBT
    [ax[0].plot(dbt_mean.iloc[i:i + 24], color='#BC204B', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[0].fill_between(np.arange(i, i + 24), dbt_min.iloc[i:i + 24], dbt_max.iloc[i:i + 24], color='#BC204B',
                        alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
    ax[0].set_ylabel(dbt_name, labelpad=2, color=tone_color)
    ax[0].yaxis.set_major_locator(MaxNLocator(7))

    # Plot DPT / RH
    [ax[1].plot(moisture_mean.iloc[i:i + 24], color='#00617F', lw=2, label='Average') for i in np.arange(0, 288)[::24]]
    [ax[1].fill_between(np.arange(i, i + 24), moisture_min.iloc[i:i + 24], moisture_max.iloc[i:i + 24], color='#00617F',
                        alpha=0.2, label='Range') for i in np.arange(0, 288)[::24]]
    ax[1].set_ylabel(moisture_name, labelpad=2, color=tone_color)
    ax[1].yaxis.set_major_locator(MaxNLocator(7))
    if not dew_point:
        ax[1].set_ylim([0, 100])

    # Plot solar
    [ax[2].plot(direct_radiation_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls='--', label=dir_rad_name) for
     i in
     np.arange(0, 288)[::24]]
    [ax[2].plot(diffuse_radiation_mean.iloc[i:i + 24], color='#FF8F1C', lw=1.5, ls=':', label=dif_rad_name) for i
     in np.arange(0, 288)[::24]]
    [ax[2].plot(global_radiation_mean.iloc[i:i + 24], color='#FF8F1C', lw=2, ls='-', label=glob_rad_name)
     for
     i in np.arange(0, 288)[::24]]
    ax[2].set_ylabel('Solar Radiation (W/m²)', labelpad=2, color=tone_color)
    ax[2].yaxis.set_major_locator(MaxNLocator(7))

    # Format plot area
    [[i.spines[spine].set_visible(False) for spine in ['top', 'right']] for i in ax]
    [[i.spines[j].set_color(tone_color) for i in ax] for j in ['bottom', 'left']]
    [i.xaxis.set_ticks(np.arange(0, 288, 24)) for i in ax]
    [i.set_xlim([0, 287]) for i in ax]
    [plt.setp(i.get_yticklabels(), color=tone_color) for i in ax]
    [i.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
                       ha='left', color=tone_color) for i in ax]
    [i.get_xaxis().set_ticklabels([]) for i in [ax[0], ax[1]]]
    [i.grid(b=True, which='major', axis='both', c=tone_color, ls='--', lw=1, alpha=0.3) for i in ax]
    [i.tick_params(length=0) for i in ax]

    ax[2].set_ylim([0, ax[2].get_ylim()[1]])

    # Legend
    handles, labels = ax[2].get_legend_handles_labels()
    lgd = ax[2].legend(bbox_to_anchor=(0.5, -0.2), loc=8, ncol=3, borderaxespad=0., frameon=False,
                       handles=[handles[0], handles[12], handles[24]], labels=[labels[0], labels[12], labels[24]])
    lgd.get_frame().set_facecolor((1, 1, 1, 0))
    [plt.setp(text, color=tone_color) for text in lgd.get_texts()]

    # Add a title
    title = plt.suptitle(
        "Monthly average diurnal profile\n{0:} - {1:} - {2:}".format(self.city, self.country, self.station_id),
        color=tone_color, y=1.025)

    # Tidy plot
    plt.tight_layout()

    # Save figure
    if save:
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "diurnal_{0:}.png".format("dew_point_temperature" if dew_point else "relative_humidity")
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
        print("Diurnal plot saved to {}".format(save_path))

    plt.close()

    return fig
