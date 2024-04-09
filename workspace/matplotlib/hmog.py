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


def plot_observable_density(ax, plot_x1s, plot_x2s, density, observations):

    ax.set_xlabel("x_1")
    ax.set_ylabel("x_2")

    # Plot for mixture model densities (panel 0)
    # Contour plot for the density evaluted at the plot_xs
    X1, X2 = np.meshgrid(plot_x1s, plot_x2s)
    density = density.reshape(X2.shape)
    ax.contourf(X1, X2, density, levels=10, cmap=density_palette)
    ax.scatter(
        observations[:, 0],
        observations[:, 1],
        color="black",
        s=5,
    )

    obs_mean = observations.mean(axis=0)
    obs_centered = observations - obs_mean

    # Step 2: Compute the covariance matrix
    cov_matrix = np.cov(obs_centered.T)

    # Step 3: Find eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # The first principal component is the eigenvector corresponding to the largest eigenvalue
    pc1 = eigenvectors[:, np.argmax(eigenvalues)]

    # Plot the first principal component as an arrow
    # Scale the arrow length for visibility; you might need to adjust the scaling factor
    arrow_start = obs_mean
    arrow_end = (
        obs_mean + pc1 * np.sqrt(np.max(eigenvalues)) * 2
    )  # Scale for visibility
    ax.quiver(
        *arrow_start,
        *(arrow_end - arrow_start),
        color="red",
        scale_units="xy",
        angles="xy",
        scale=1,
    )


def plot_cross_entropy(ax, cross_entropy):
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Cross Entropy")
    ax.plot(cross_entropy, color="black")


def plot_latent_densities(ax, plot_ys, true_latent_density, latent_density):
    ax.set_xlabel("y")
    ax.set_ylabel("Density")
    # Plot for mixture model densities (panel 0)
    ax.plot(plot_ys, true_latent_density, label="True Latent Density", color="black")
    ax.plot(plot_ys, latent_density, label="Learned Latent Density", color="red")


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("hmog.json"), "r") as file:
        data = json.load(file)

    fig = plt.figure(figsize=(12, 8))
    outer_grid = gridspec.GridSpec(2, 1, height_ratios=[2, 1])

    # Top Row (A, B, C)
    top_row = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_grid[0])
    axA = plt.subplot(top_row[0])
    axB = plt.subplot(top_row[1])
    axC = plt.subplot(top_row[2])

    # Right Column (B and D)
    bottom_row = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_grid[1])

    axD = plt.subplot(bottom_row[0])
    axE = plt.subplot(bottom_row[1])
    axF = plt.subplot(bottom_row[2])

    for i, ax in enumerate([axA, axB, axC]):
        plot_x1s = data["plot-range-x1"]
        plot_x2s = data["plot-range-x2"]
        density = np.array(data["observable-densities"][i])
        observations = np.array(data["observations"])
        plot_observable_density(ax, plot_x1s, plot_x2s, density, observations)

    cross_entropy = -np.array(data["log-likelihoods"])

    plot_cross_entropy(axD, cross_entropy)

    true_latent_density = np.array(data["mixture-densities"][0])

    for i, ax in enumerate([axE, axF]):
        plot_ys = data["plot-range-y"]
        latent_density = np.array(data["mixture-densities"][i + 1])
        plot_latent_densities(ax, plot_ys, true_latent_density, latent_density)

    plt.savefig(get_plot_path("hmog.png"))
