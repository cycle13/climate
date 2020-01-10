from climate.common.constants import time_of_year_mask
from climate.common.helpers import angle_between

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def generate_radiation_rose_values(sky_matrix, patch_vectors, season_period="Annual", day_period="Daily", nsector=36):
    # Filter and summate passed sky matrix by time (removing hours not to be included)
    mask = np.array([time_of_year_mask[day_period], time_of_year_mask[season_period]]).all(axis=0)
    filtered_rad = sky_matrix[mask]
    cumulative_rad = filtered_rad.sum(axis=0)

    # Calculate vectors into which patches will be binned
    rose_angles = np.radians(np.linspace(0, 360, nsector, endpoint=False))
    radiation_rose_vectors = np.stack([-np.sin(rose_angles), np.cos(rose_angles), np.zeros(len(rose_angles))]).T

    # Get the angle between each patch centroid and rose vector
    vector_angles_top = []
    for radiation_rose_vector in radiation_rose_vectors:
        vector_angles_bottom = []
        for patch_vector in patch_vectors:
            vector_angles_bottom.append(angle_between(patch_vector, radiation_rose_vector, degrees=False))
        vector_angles_top.append(vector_angles_bottom)
    patch_rose_vector_angles = np.array(vector_angles_top)

    # Remove radiation beyond view of the radiation rose vector
    te = cumulative_rad * np.cos(patch_rose_vector_angles)
    radiation_rose_values = np.where(te.sum(axis=1) < 0, 0, te.sum(axis=1))
    radiation_rose_values = np.roll(radiation_rose_values, -1)

    return rose_angles, radiation_rose_values


def radiation_rose(self, season_period="Annual", day_period="Daily", n_sector=36, cmap=None, tone_color="k", same_scale=False, save=False):
    # Create values to plot
    rose_angles, dir_rad = generate_radiation_rose_values(self.direct_sky_matrix, self.patch_vectors, season_period=season_period,
                                                          day_period=day_period, nsector=n_sector)
    rose_angles, dif_rad = generate_radiation_rose_values(self.diffuse_sky_matrix, self.patch_vectors, season_period=season_period,
                                                          day_period=day_period, nsector=n_sector)
    tot_rad = dir_rad + dif_rad

    max_val = tot_rad.max() / 1000

    figs = []
    for tx, vals in {"Direct Radiation": dir_rad,
                     "Diffuse Radiation": dif_rad,
                     "Total Radiation": tot_rad}.items():

        # Construct the save_path and create directory if it doesn't exist
        save_path = self.file_path.parent / "{}_Plot".format(self.file_path.stem) / "radiationrose_{}_{}{}.png".format(
            tx.replace(" Radiation", ""), season_period, day_period)

        vals /= 1000
        colors = [cm.get_cmap('OrRd')(i) if cmap is None else cm.get_cmap(cmap)(i) for i in np.interp(vals, [min(vals), max(vals)], [0, 1])]
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), subplot_kw={'projection': "polar"})
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.bar(
            rose_angles[::-1],
            vals,
            width=(np.pi * 2 / (n_sector + 1)) * 0.95, zorder=5, bottom=0.0, color=colors, alpha=1, edgecolor="w",
            linewidth=0.1)
        if same_scale:
            ax.set_ylim(0, max_val)
        ax.spines['polar'].set_visible(False)
        ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"], color=tone_color)
        plt.setp(ax.get_yticklabels(), color=tone_color)
        ax.set_title("{} (W/m²) - {} - {}\n{} - {} - {}".format(tx, season_period, day_period, self.city, self.country, self.station_id), color=tone_color, loc="center", va="bottom", ha="center", fontsize="large", y=1)
        plt.tight_layout()
        figs.append(fig)

        if save:
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, bbox_inches="tight", dpi=300, transparent=False)
            print("{0:} rose saved to {1:}".format(tx, save_path))

        plt.close()

    return figs
