### Imports ###
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import seaborn as sns

# Local

from common import get_result_path, get_plot_path

### Globals ###

density_palette = sns.color_palette("twilight", as_cmap=True)

density_range = 0.2

### Plotting Functions ###


def plot_simplex(ax, plot_xs, plot_ys, vertices, density, true_cat, avg_obs=None):
    # Original simplex coordinates
    # Contour plot for the density evaluted at the plot_xs
    X, Y = np.meshgrid(plot_xs, plot_ys)
    density = np.array(density).reshape(Y.shape)
    heatmap = ax.contour(
        X, Y, density, levels=10, cmap=density_palette, vmin=0, vmax=0.2
    )

    ax.set_aspect("equal")
    ax.axis("off")
    # Draw the simplex boundary (equilateral triangle), adjusted for centering
    triangle = Polygon(vertices, edgecolor="k", fill=False, linewidth=1.5, zorder=10)
    ax.add_patch(triangle)

    ax.scatter(
        true_cat[0], true_cat[1], c="black", s=50, label="True Categorical", zorder=20
    )
    if avg_obs is not None:
        ax.scatter(
            avg_obs[0],
            avg_obs[1],
            c="red",
            s=100,
            label="Average Observations",
            zorder=30,
        )

    return heatmap


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("categorical-inference.json"), "r") as file:
        data = json.load(file)

    naxs = 5
    fig, axs = plt.subplots(1, naxs, figsize=(20, 5), width_ratios=[1, 1, 1, 1, 0.1])

    true_categorical = data["true-categorical"]
    observations = data["observations"]
    average_observations = [None] + data["average-observations"]
    dirichlets = data["dirichlets"]
    plot_xs = data["plot-xs"]
    plot_ys = data["plot-ys"]
    vertices = data["vertices"]
    dirichlet_densities = data["dirichlet-densities"]

    # all but the last plot

    for i, ax in enumerate(axs[:-1]):
        heatmap = plot_simplex(
            ax,
            plot_xs,
            plot_ys,
            vertices,
            dirichlet_densities[i],
            true_categorical,
            average_observations[i],
        )
        if i == naxs - 2:
            fax = axs[naxs - 1]
            # stretch color bar across panel
            fig.colorbar(heatmap, ax=ax, cax=fax)

    plt.tight_layout()
    plt.savefig(get_plot_path("categorical-inference.png"))
