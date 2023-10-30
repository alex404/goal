import matplotlib.pyplot as plt
import json

# Local
from common import get_result_path, get_plot_path

# Define the filenames based on the distributions
filenames = [
    "univariate-binomial.json",
    "univariate-categorical.json",
    "univariate-poisson.json",
    "univariate-normal.json",
]

# Create a 2x2 subplot layout
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
axs = axs.ravel()

for idx, fname in enumerate(filenames):
    full_path = get_result_path(fname)

    # Read the JSON file
    with open(full_path, "r") as f:
        data = json.load(f)

    # Extract and process the data
    true_density = data["true-density"]
    source_density = data["source-density"]
    natural_density = data["natural-density"]
    samples = data["samples"]
    range_data = data["range"]

        # Flatten samples for univariate-normal data
    if fname == "univariate-normal.json":
        samples = [item for sublist in samples for item in sublist]

    # Plot the data on the respective subplot
    ax = axs[idx]

    # Plot histogram (using numpy to compute the histogram)
    n, bins, patches = ax.hist(
        samples, bins="auto", color="black", label="Samples", alpha=0.5, density=True
    )

    # If samples are integers, adjust the x-axis ticks
    if all(isinstance(sample, int) for sample in samples):
        ax.set_xticks(bins)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Plot densities
    ax.plot(range_data, true_density, color="black", lw=2, label="True Density")
    ax.plot(range_data, source_density, color="blue", lw=2, label="Source Density")
    ax.plot(
        range_data,
        natural_density,
        color="red",
        lw=2,
        linestyle="dashed",
        label="Natural Density",
    )

    # Labels and title
    ax.set_xlabel("X")
    ax.set_ylabel("Density")
    ax.set_title(data["title"])

    # Setting legend on the last plot for space efficiency
    if idx == 3:
        ax.legend()

plt.tight_layout()

plot_file_path = get_plot_path("univariate-analytical.png")

# Save the plot to a file
plt.savefig(plot_file_path)
