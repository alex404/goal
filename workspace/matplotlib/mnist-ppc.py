import json
import numpy as np
import matplotlib.pyplot as plt
from common import get_result_path,get_plot_path  # Ensure this is correctly pointing to your utility functions

# Load the Gaussian-Boltzmann data
filename = "mnist-ppc.json"  # Update this to your actual filename
with open(get_result_path(filename), "r") as file:
    gb_data = json.load(file)

# Extract the data
xys = np.array(gb_data["density-xys"])
initial_density = np.array(gb_data["initial-density"])
learned_density = np.array(gb_data["learned-density"])

# Determine the number of samples per side of the grid
num_samples_side = int(np.sqrt(len(xys)))

# Since xys is a list of tuples (flattened grid), separate it into two arrays for X and Y coordinates
X = xys[:, 0].reshape((num_samples_side, num_samples_side))
Y = xys[:, 1].reshape((num_samples_side, num_samples_side))
initial_density = initial_density.reshape((num_samples_side, num_samples_side))
learned_density = learned_density.reshape((num_samples_side, num_samples_side))

# Plotting
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), sharex=True, sharey=True)

# Initial density contour plot
ctr_initial = axes[0].contourf(X, Y, initial_density, levels=20, cmap='viridis')
fig.colorbar(ctr_initial, ax=axes[0])
axes[0].set_title('Initial Density')
axes[0].set_xlabel('X Coordinate')
axes[0].set_ylabel('Y Coordinate')

# Learned density contour plot
ctr_learned = axes[1].contourf(X, Y, learned_density, levels=20, cmap='viridis')
fig.colorbar(ctr_learned, ax=axes[1])
axes[1].set_title('Learned Density')
axes[1].set_xlabel('X Coordinate')
axes[1].set_ylabel('Y Coordinate')

plt.tight_layout()
plot_file_path = get_plot_path("mnist-ppc.png")
plt.savefig(plot_file_path)
