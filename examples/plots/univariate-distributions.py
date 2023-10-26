import matplotlib.pyplot as plt
import json
import os

def plot_from_json(directory):
    # Define the filenames based on the distributions
    filenames = ["univariate-binomial.json", "univariate-categorical.json", "univariate-poisson.json", "univariate-normal.json"]

    # Create a 2x2 subplot layout
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    axs = axs.ravel()

    for idx, fname in enumerate(filenames):
        full_path = os.path.join(directory, fname)

        # Read the JSON file
        with open(full_path, 'r') as f:
            data = json.load(f)

        # Extract and process the data
        histogram_bins = data["histogram-bins"]
        histogram_weights = data["histogram-weights"]
        
        true_density = data["true-density"]
        source_density = data["source-density"]
        natural_density = data["natural-density"]

        range_data = data["range"]
        
        # Plot the data on the respective subplot
        ax = axs[idx]

        # Plot histogram
        ax.bar(histogram_bins, histogram_weights, width=0.8, color='black', label="Samples", alpha=0.5)

        # Plot densities
        ax.plot(range_data, true_density, color='black', lw=2, label="True Density")
        ax.plot(range_data, source_density, color='blue', lw=2, label="Source Density")
        ax.plot(range_data, natural_density, color='red', lw=2, linestyle='dashed', label="Natural Density")

        # Labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Density')
        ax.set_title(data["title"])
        
        # Setting legend on the last plot for space efficiency
        if idx == 3:
            ax.legend()

    plt.tight_layout()
    output_file = os.path.join(directory, "all_distributions.png")
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    import sys
    plot_from_json(sys.argv[1])

