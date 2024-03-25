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
    confidence_ellipses = data["component-confidence-ellipses"]
    for i, ellipse in enumerate(confidence_ellipses):
        xs, ys = zip(*ellipse)
        ax.plot(xs, ys, label=f"Component {i}")


def plot_observable_density(ax, data):
    ax.set_xlabel("x_1")
    ax.set_ylabel("x_2")

    # Plot for mixture model densities (panel 0)
    plot_xs = data["plot-range-x"]
    plot_ys = data["plot-range-y"]
    density = np.array(data["learned-density"])
    obss = np.array(data["posterior-observations"])
    # Contour plot for the density evaluted at the plot_xs
    X1, X2 = np.meshgrid(plot_xs, plot_ys)
    density = density.reshape(X2.shape)
    ax.contourf(X1, X2, density, levels=100, cmap="viridis")
    for pclr, (obs1, obs2) in zip(posterior_colours, obss):
        ax.scatter(obs1, obs2, color=pclr, s=50, label="Observations")


def plot_prior(ax, data):

    second_moment_matrix = np.array(data["prior-moment-matrix"][0])
    heatmap = ax.imshow(
        second_moment_matrix, cmap="viridis", interpolation="nearest", vmin=0, vmax=1
    )
    ax.set_xlabel("Neuron ID")
    ax.set_ylabel("Neuron ID")
    ax.set_xticks(np.arange(second_moment_matrix.shape[1]))
    ax.set_yticks(np.arange(second_moment_matrix.shape[0]))
    return heatmap


def plot_posterior(ax, data):
    posterior_matrices = np.array(data["posterior-moment-matrices"][0])

    # Create a nested GridSpec for 2x2 subplots within axD
    inner_grid = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=ax, wspace=0.1, hspace=0.1
    )
    heatmaps = []

    for i, matrix in enumerate(posterior_matrices):
        # Calculate the position of the current subplot
        inner_ax = plt.subplot(inner_grid[i])

        # Plot the heatmap
        heatmap = inner_ax.imshow(
            matrix, cmap="viridis", interpolation="nearest", vmin=0, vmax=1
        )
        heatmaps.append(heatmap)
        inner_ax.set_xticks([])
        inner_ax.set_yticks([])
        rect = Rectangle(
            (0, 0),
            1,
            1,
            linewidth=5,
            edgecolor=posterior_colours[i],
            facecolor="none",
            transform=inner_ax.transAxes,
        )
        inner_ax.add_patch(rect)

    ax.axis("off")
    return heatmaps


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("gaussian-boltzmann.json"), "r") as file:
        data = json.load(file)

    fig = plt.figure(figsize=(12, 10))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.2, width_ratios=[1, 1])

    # Left Column (A and C)
    left_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[0], hspace=0.1
    )
    axA = plt.subplot(left_column[0])
    axC = plt.subplot(left_column[1], sharex=axA, sharey=axA)

    # Right Column (B and D)
    right_column = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer_grid[1], hspace=0.1
    )
    axB = plt.subplot(right_column[0])
    axD = plt.subplot(right_column[1], sharex=axB, sharey=axB)

    # Example plotting calls
    # Replace these with your plot functions
    plot_likelihood(axA, data)
    heatmap = plot_prior(axB, data)
    fig.colorbar(heatmap, ax=axB)
    plot_observable_density(axC, data)
    plot_posterior(axD, data)

    # Adjustments

    plt.savefig(get_plot_path("gaussian-boltzmann.png"))
