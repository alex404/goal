### Imports ###

import matplotlib.pyplot as plt
import json
import seaborn as sns

# Local
from common import get_result_path, get_plot_path

### Globals ###

nbins = 20

### Plotting Functions ###


def plot_von_mises(ax, data):

    # Extract and process the data
    samples = data["samples"]
    num_samples = len(samples)

    # Define bins
    pi = 3.141592653589793
    bin_width = (2 * pi) / nbins  # 2π divided by 20
    bins = [i * bin_width for i in range(nbins + 1)]  # 21 edges for 20 bins
    ax.set_xlim([0, 2 * pi])  # 2π
    ax.set_xticks([i * pi / 2 for i in range(5)])  # Ticks for 0, π/4, π/2, ..., 2π
    ax.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
    bar_width = ((2 * pi) / nbins) / (num_samples + 0.5)

    plot_general(ax, data, bins, bar_width)


def plot_com_poisson(ax, data):

    samples = data["samples"]
    num_samples = len(samples)
    range_data = data["range"]
    bins = range(min(range_data), max(range_data) + 2)  # Creating integer bins
    bar_width = 1.0 / (num_samples + 0.5)

    plot_general(ax, data, bins, bar_width)


def plot_general(ax, data, bins, bar_width):
    # Plot multiple sample histograms
    title = data["title"]
    true_density = data["true-density"]
    natural_density = data["natural-density"]
    samples = data["samples"]
    range_data = data["range"]
    num_samples = len(samples)
    colors = sns.color_palette("Set2", num_samples)

    for dataset_idx, sample_set in enumerate(samples):
        # Adjust position for each histogram dataset
        adjusted_bins = [x + 0.025 + dataset_idx * bar_width for x in bins]

        # Plot histogram with adjusted position
        ax.hist(
            sample_set,
            bins=adjusted_bins,
            color=colors[dataset_idx],
            alpha=0.9,
            density=True,
            label=f"Sample Size {len(sample_set)}",
            width=bar_width,
        )

    # Plot densities
    ax.plot(range_data, true_density, color="black", lw=3, label="True Density")
    ax.plot(
        range_data,
        natural_density,
        color="red",
        lw=3,
        linestyle="dashed",
        label="Fit Density",
    )

    # Labels and title
    ax.set_xlabel("x")
    ax.set_ylabel("Density")
    ax.xaxis.grid(True)
    ax.legend()
    # Caps first letter and bold title
    ax.set_title(title.capitalize(), fontweight="bold")


def plot_gamma(ax, data):

    samples = data["samples"]
    num_samples = len(samples)
    range_data = data["range"]
    bin_width = (max(range_data) - min(range_data)) / nbins
    bar_width = 1.0 / (num_samples + 0.5)
    bins = [i * bin_width for i in range(nbins + 1)]  # 21 edges for 20 bins

    plot_general(ax, data, bins, bar_width)


### Main ###

if __name__ == "__main__":

    # Load the data
    with open(get_result_path("univariate-von-mises.json"), "r") as file:
        von_mises_data = json.load(file)
    with open(get_result_path("univariate-com-poisson.json"), "r") as file:
        com_poisson_data = json.load(file)
    with open(get_result_path("univariate-gamma.json"), "r") as file:
        gamma_data = json.load(file)

    fig, axes = plt.subplots(3, 1, figsize=(8, 12), constrained_layout=True)

    # Example plotting calls
    plot_von_mises(axes[0], von_mises_data)
    plot_com_poisson(axes[1], com_poisson_data)
    plot_gamma(axes[2], gamma_data)

    plt.savefig(get_plot_path("univariate-numerical.png"))
