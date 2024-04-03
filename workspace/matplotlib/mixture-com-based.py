### Imports ###

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

# Local
from common import get_result_path, get_plot_path


### Globals ###

contour_resolution = 100
scatter_palette = sns.color_palette("Set1", n_colors=2)
correlation_palette = sns.color_palette("RdBu", as_cmap=True)


### Functions ###


## Helper Functions ##


## Plotting Functions ##


def plot_mean_scatter(ax, data):

    smp_ffs = np.array(data["sample-fano-factors"])
    cbm_ffs = np.array(data["cbm-fano-factors"])
    data_min = min(smp_ffs.min(), cbm_ffs.min())
    data_max = max(smp_ffs.max(), cbm_ffs.max())
    # Scatter plot
    padding = (data_max - data_min) * 0.05
    data_range = [data_min - padding, data_max + padding]
    ax.scatter(smp_ffs, cbm_ffs, s=10, color=scatter_palette[0])
    ax.set_xlim(data_range)
    ax.set_ylim(data_range)
    ax.set_aspect("equal", "box")
    ax.set_ylabel("Sample Means")
    ax.set_xlabel("CoM-Based Mixture Means")


def plot_covariance_scatter(ax, data):

    smp_cvrs = np.array(data["sample-covariances"])
    cbm_cvrs = np.array(data["cbm-covariances"])
    data_min = min(smp_cvrs.min(), cbm_cvrs.min())
    data_max = max(smp_cvrs.max(), cbm_cvrs.max())
    # Scatter plot
    padding = (data_max - data_min) * 0.05
    data_range = [data_min - padding, data_max + padding]
    ax.scatter(smp_cvrs, cbm_cvrs, s=10, color=scatter_palette[1])
    ax.set_xlim(data_range)
    ax.set_ylim(data_range)
    ax.set_aspect("equal", "box")
    ax.set_ylabel("Sample Covariances")
    ax.set_xlabel("CoM-Based Mixture Covariances")


def plot_cbm_correlations(ax, data, correlation_palette="seismic"):
    correlations = np.array(data["cbm-correlation-matrix"])

    # Perform hierarchical clustering to find a new ordering for the matrix
    linkage_matrix = linkage(correlations, method="ward")
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = leaves_list(linkage_matrix)

    # Reorder the correlation matrix according to the clustering result
    sorted_correlations = correlations[order, :][:, order]

    # Plotting the sorted correlation matrix
    heatmap = ax.imshow(
        sorted_correlations,
        cmap=correlation_palette,
        interpolation="nearest",
        vmin=-1,
        vmax=1,
    )

    ax.set_xlabel("Observable Dimension")
    ax.set_ylabel("Observable Dimension")
    ax.set_xticks(np.arange(sorted_correlations.shape[1]))
    ax.set_yticks(np.arange(sorted_correlations.shape[0]))
    # Optional: Set tick labels to the sorted order (or any meaningful labels)
    ax.set_xticklabels(order)
    ax.set_yticklabels(order)

    return heatmap


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("mixture-com-based.json"), "r") as file:
        data = json.load(file)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    # Example plotting calls
    # Replace these with your plot functions
    plot_mean_scatter(axes[0], data)
    plot_covariance_scatter(axes[1], data)
    # heatmap = plot_correlations(axes[1], data)
    # fig.colorbar(heatmap, ax=axes[1])
    heatmap = plot_cbm_correlations(axes[2], data)
    fig.colorbar(heatmap, ax=axes[2])

    # Adjustments

    plt.savefig(get_plot_path("mixture-com-based.png"))
