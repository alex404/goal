### Imports ###

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import matplotlib.gridspec as gridspec

# Local
from common import get_result_path, get_plot_path


### Globals ###

contour_resolution = 100
posterior_colours = sns.color_palette("Set1", n_colors=4)
density_palette = sns.color_palette("twilight", as_cmap=True)


### Functions ###


## Helper Functions ##


## Plotting Functions ##


def plot_observable_density(ax, data):
    ax.set_xlabel("x_1")
    ax.set_ylabel("x_2")

    # Plot for mixture model densities (panel 0)
    plot_x1s = data["plot-range-x1"]
    plot_x2s = data["plot-range-x2"]
    density = np.array(data["true-observable-density"])
    # Contour plot for the density evaluted at the plot_xs
    X1, X2 = np.meshgrid(plot_x1s, plot_x2s)
    density = density.reshape(X2.shape)
    ax.contour(X1, X2, density, levels=10, colors="black")


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("hmog.json"), "r") as file:
        data = json.load(file)

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))

    # # Right Column (B and D)
    # right_column = gridspec.GridSpecFromSubplotSpec(
    #     2, 1, subplot_spec=outer_grid[1], hspace=0.1
    # )
    # axB = plt.subplot(right_column[0])
    # axD = plt.subplot(right_column[1], sharex=axB, sharey=axB)
    #
    # # Example plotting calls
    # # Replace these with your plot functions
    # plot_likelihood(axA, data)
    # heatmap = plot_prior(axB, data)
    # fig.colorbar(heatmap, ax=axB)
    plot_observable_density(ax, data)
    # plot_posterior(axD, data)
    #
    # Adjustments

    plt.savefig(get_plot_path("hmog.png"))
