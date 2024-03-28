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


### Functions ###


## Helper Functions ##


## Plotting Functions ##


def plot_likelihood(ax, data):
    # Set global configurations for both plots
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("Orientation")
    ax.set_ylabel("Rate")
    ax.set_ylim(0, 8)
    ax2 = ax.twinx()
    ax2.set_ylabel("Total Rate")
    ax2.set_ylim(0, 30)
    ax2.spines["right"].set_visible(True)

    # Plot for affine population
    ax2.plot(data["plot-zs"], data["sum-of-tuning-curves"], lw=4)
    ax2.plot(data["plot-zs"], data["regression"], lw=4, linestyle="--")
    for curve in data["tuning-curves"]:
        ax.plot(data["plot-zs"], curve, lw=4, color="black")


def plot_observable_density(ax, data):
    correlation_matrix = np.array(data["correlation-matrix"])
    heatmap = ax.imshow(
        correlation_matrix, cmap="viridis", interpolation="nearest", vmin=-1, vmax=1
    )
    ax.set_xlabel("Neuron ID")
    ax.set_ylabel("Neuron ID")
    ax.set_xticks(np.arange(correlation_matrix.shape[1]))
    ax.set_yticks(np.arange(correlation_matrix.shape[0]))
    return heatmap


def plot_prior(ax, data):
    ax.plot(data["plot-zs"], data["prior-density"], lw=4)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("Orientation")
    ax.set_ylabel("Density")
    return


def plot_posterior(ax, data):
    return


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("population-code-von-mises.json"), "r") as file:
        data = json.load(file)

    fig = plt.figure(figsize=(12, 8))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.4, hspace=0.2, width_ratios=[1, 1])

    # Left Column (A and C)
    left_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[0], hspace=0.1
    )
    axA = plt.subplot(left_column[0])
    axC = plt.subplot(left_column[1])

    # Right Column (B and D)
    right_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[1], hspace=0.1
    )
    axB = plt.subplot(right_column[0])
    axD = plt.subplot(right_column[1], sharex=axB, sharey=axB)

    # Example plotting calls
    # Replace these with your plot functions
    plot_likelihood(axA, data)
    plot_prior(axB, data)
    heatmap = plot_observable_density(axC, data)
    fig.colorbar(heatmap, ax=axC, label="Correlation")
    plot_posterior(axD, data)

    # Adjustments

    plt.savefig(get_plot_path("population-code-von-mises.png"))
