import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from climate.common.constants import time_of_year_mask
from climate.common.helpers import angle_between


def plot_radiation_rose(self, save=False):

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
