### Imports ###

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sns
import matplotlib.gridspec as gridspec

# Local
from common import get_result_path, get_plot_path


### Globals ###

simplex_resolution = 500
vertices = np.array([[0, 0], [1, 0], [0.5, 1]])
colors = sns.color_palette("Set1", n_colors=3)


### Functions ###


## Helper Functions ##


def barycentric_to_cartesian(barycentric):
    """Converts a point from barycentric to Cartesian coordinates."""
    return np.dot(vertices.T, barycentric)


## Plotting Functions ##


def plot_likelihood(ax, data):
    ax.set_xlabel("x")
    ax.set_ylabel("Density")
    # Plot for mixture model densities (panel 0)
    plot_samples = data["plot-samples"]

    # Plot the individual component densities
    for comp_density, color in zip(
        data["component-densities"],
        colors,
    ):
        ax.plot(plot_samples, comp_density, label=f"Component Density", color=color)


def plot_observable_density(ax, data):
    ax.set_xlabel("x")
    ax.set_ylabel("Density")

    # Plot for mixture model densities (panel 0)
    plot_samples = data["plot-samples"]
    mix_densities = data["mixture-density"]
    ax.plot(plot_samples, mix_densities, label="Mixture Density", color="black")


def plot_prior(ax, data):
    mixture_weights = data["weights"]

    # Create a mesh grid, adjusted to center the triangle
    x = np.linspace(0, 1, simplex_resolution)
    y = np.linspace(0, 1, simplex_resolution)
    X, Y = np.meshgrid(x, y)

    # Mask for the simplex area (equilateral triangle), adjusted for centering
    mask = (Y <= 2 * X) & (Y <= ((1 - 2 * (X - 0.5))))

    # Initialize an RGBA image
    image = np.zeros((simplex_resolution, simplex_resolution, 4))

    # Calculate RGB values based on barycentric coordinates
    for i in range(simplex_resolution):
        for j in range(simplex_resolution):
            if mask[j, i]:
                # Adjust calculations to account for the centering of the triangle
                barycentric = np.array([1 - X[j, i] - Y[j, i], X[j, i], Y[j, i]])
                color = np.dot(np.array([colors[k] for k in range(3)]).T, barycentric)
                image[j, i] = np.append(color, 1)

    # Display the image, correctly orienting and scaling it
    ax.imshow(image, extent=(0, 1, 0, 1), origin="lower", aspect="auto")

    # Draw the simplex boundary (equilateral triangle), adjusted for centering
    triangle = Polygon(vertices, edgecolor="k", fill=False, linewidth=1.5, zorder=10)
    ax.add_patch(triangle)

    x, y = barycentric_to_cartesian(mixture_weights)
    # Add large black point for the mixture weights in barycentric coordinates
    ax.plot(x, y, "wo", markersize=5, zorder=20)
    ax.plot(x, y, "ko", markersize=3, zorder=30)

    # Configure plot to ensure the triangle is centered and fully visible
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)


def plot_posterior(ax, data):
    posteriors = data["posteriors"]

    # Draw the simplex boundary (equilateral triangle), adjusted for centering
    triangle = Polygon(vertices, edgecolor="k", fill=False, linewidth=1.5, zorder=10)
    ax.add_patch(triangle)

    # Plot each posterior weight as points
    for posterior, color in zip(posteriors, colors):
        x, y = barycentric_to_cartesian(posterior)
        ax.plot(
            x, y, "o", c=color, markersize=8, zorder=30, alpha=0.8
        )  # Adjust alpha for visibility

    # Configuration
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("mixture-normal.json"), "r") as file:
        data = json.load(file)

    fig = plt.figure(figsize=(12, 6))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.2, width_ratios=[1.75, 1])

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
    plot_prior(axB, data)
    plot_observable_density(axC, data)
    plot_posterior(axD, data)

    # Adjustments
    axA.tick_params(labelbottom=False)
    axB.tick_params(labelbottom=False)

    plt.savefig(get_plot_path("mixture-normal.png"))
