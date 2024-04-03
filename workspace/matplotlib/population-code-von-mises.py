### Imports ###

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

# Local
from common import get_result_path, get_plot_path


### Globals ###

contour_resolution = 100
posterior_colours = sns.color_palette("Set1", n_colors=4)
correlation_palette = sns.color_palette("seismic_r", as_cmap=True)

### Functions ###


## Helper Functions ##


## Plotting Functions ##


def plot_likelihood(ax, data):
    # Set global configurations for both plots
    ax2 = ax.twinx()
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("Orientation")
    ax.set_ylabel("Rate")
    ax.set_ylim(0, 4)
    ax.set_yticks([0, 2, 4])
    ax2.set_ylabel("Total Rate")
    ax2.set_ylim(0, 10)
    ax2.set_yticks([0, 5, 10])
    ax2.spines["right"].set_visible(True)

    # Plot for affine population
    ax2.plot(data["plot-zs"], data["sum-of-tuning-curves"], lw=4)
    ax2.plot(data["plot-zs"], data["regression"], lw=4, linestyle="--")
    for curve in data["tuning-curves"]:
        ax.plot(data["plot-zs"], curve, lw=4, color="black")


def plot_observable_density(ax, data):
    correlation_matrix = np.array(data["correlation-matrix"])
    heatmap = ax.imshow(
        correlation_matrix,
        cmap=correlation_palette,
        interpolation="nearest",
        vmin=-1,
        vmax=1,
    )
    # Number of neurons (assumed to be even and equal to the dimension of the correlation matrix)
    n_neurons = correlation_matrix.shape[0]
    div = n_neurons // 2

    # Generating preferred stimuli labels uniformly sampled over the circle
    # The labels are represented as fractions of π
    labels = [f"${{\\frac{{{i}\\pi}}{{{div}}}}}$" for i in range(0, n_neurons)]

    # Setting the tick labels for both axes
    ax.set_xticks(np.arange(n_neurons))
    ax.set_yticks(np.arange(n_neurons))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.spines["right"].set_visible(True)
    ax.spines["top"].set_visible(True)

    # Setting axis labels
    ax.set_xlabel("Preferred Stimulus")
    ax.set_ylabel("Preferred Stimulus")

    return heatmap


def plot_prior(ax, data):
    ax.plot(data["plot-zs"], data["prior-density"], lw=4)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("Orientation")
    ax.set_ylabel("Density")


def plot_posterior(ax, data):
    posteriors = np.array(data["posterior-densities"])
    for i, posterior in enumerate(posteriors):
        ax.plot(data["plot-zs"], posterior, lw=4, color=posterior_colours[i])
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel("Orientation")
    ax.set_ylabel("Density")

    ax2 = ax.twinx()
    ax2.set_ylabel("Spike Count")
    ax2.set_ylim(0, 10)
    ax2.set_yticks([0, 5, 10])
    ax2.spines["right"].set_visible(True)

    # Plot for affine population
    stims = data["preferred-stimuli"] + [2 * np.pi]
    for i, observation in enumerate(data["observations"]):
        observation = observation + [observation[0]]
        ax2.plot(stims, observation, lw=2, color=posterior_colours[i], dashes=[2, 2])
        ax2.scatter(
            stims,
            observation,
            s=100,
            color=posterior_colours[i],
            marker="o",
            facecolors="white",
            lw=2,
            zorder=10,
        )


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("population-code-von-mises.json"), "r") as file:
        data = json.load(file)

    fig = plt.figure(figsize=(12, 8))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.3, width_ratios=[1, 1])

    # Left Column (A and C)
    left_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[0], hspace=0.25
    )
    axA = plt.subplot(left_column[0])
    axC = plt.subplot(left_column[1])

    # Right Column (B and D)
    right_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[1], hspace=0.25
    )
    axB = plt.subplot(right_column[0])
    axD = plt.subplot(right_column[1], sharex=axB, sharey=axB)

    # Example plotting calls
    # Replace these with your plot functions
    plot_likelihood(axA, data)
    plot_prior(axB, data)
    heatmap = plot_observable_density(axC, data)
    fig.colorbar(heatmap, ax=axC, label="Correlation", ticks=[-1, 0, 1])
    plot_posterior(axD, data)

    # Adjustments

    plt.savefig(get_plot_path("population-code-von-mises.png"))
