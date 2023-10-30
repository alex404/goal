import matplotlib.pyplot as plt
import json
import seaborn as sns

# Local
from common import get_result_path, get_plot_path

# Define the filenames based on the distributions
filenames = [
    "univariate-von-mises.json",
    "univariate-com-poisson.json"
]

# Analyze the number of sample datasets in the first file to set the color palette
with open(get_result_path(filenames[0]), "r") as f:
    data = json.load(f)
num_samples = len(data["samples"])

colors = sns.color_palette("Set2", num_samples)

# Create a 1x2 subplot layout for the two distributions
fig, axs = plt.subplots(2, 1, figsize=(8, 8))
axs = axs.ravel()

for idx, fname in enumerate(filenames):
    full_path = get_result_path(fname)

    # Read the JSON file
    with open(full_path, "r") as f:
        data = json.load(f)

    # Extract and process the data
    true_density = data["true-density"]
    natural_density = data["natural-density"]
    samples = data["samples"]
    range_data = data["range"]

    # Plot the data on the respective subplot
    ax = axs[idx]

    # Define bins
    if "von-mises" in fname:
        nbins = 20
        pi = 3.141592653589793
        bin_width = (2 * pi) / nbins  # 2π divided by 20
        bins = [i * bin_width for i in range(nbins + 1)]  # 21 edges for 20 bins
        ax.set_xlim([0, 2 * pi]) # 2π
        ax.set_xticks([i * pi / 2 for i in range(5)]) # Ticks for 0, π/4, π/2, ..., 2π
        ax.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
        bar_width = ((2*pi)/nbins) / (num_samples + 0.5)

    else:
        bins = range(min(range_data), max(range_data) + 2)  # Creating integer bins
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        bar_width = 1.0 / (num_samples + 0.5)

    # Plot multiple sample histograms
    for dataset_idx, sample_set in enumerate(samples):
        # Adjust position for each histogram dataset
        adjusted_bins = [x + 0.025 + dataset_idx * bar_width for x in bins]

        # Plot histogram with adjusted position
        ax.hist(sample_set, bins=adjusted_bins, color=colors[dataset_idx], alpha=0.9, density=True, 
                label=f"Samples {len(sample_set)}", width=bar_width)

    # Plot densities
    ax.plot(range_data, true_density, color="black", lw=3, label="True Density")
    ax.plot( range_data, natural_density, color="red", lw=3, linestyle="dashed", label="Fit Density")

    # Labels and title
    ax.set_xlabel("X")
    ax.set_ylabel("Density")
    ax.set_title(data["title"])
    ax.xaxis.grid(True)

    # Setting legend on the second plot for space efficiency
    if idx == 0:
        ax.legend()

plt.tight_layout()

plot_file_path = get_plot_path("univariate-numerical.png")

# Save the plot to a file
plt.savefig(plot_file_path)

