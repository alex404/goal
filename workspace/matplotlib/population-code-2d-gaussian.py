import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# Local
from common import get_result_path, get_plot_path

### Load and prepare the data
with open(get_result_path("population-code-2d-gaussian.json"), "r") as file:
    data = json.load(file)

# Extracting the necessary elements from the data
xys = np.array(data["xys"])
stcs = np.array(data["sum-of-tuning-curves"])
stcerrs = np.array(data["estimation-difference"])
pss = np.array(data["preferred-stimuli"])
true_density = np.array(data["true-density"])
initial_density = np.array(data["initial-density"])
learned_density = np.array(data["learned-density"])
regmn, regmx = data["regression-bounds"]

# Prepare the grid
unique_xs = np.unique(xys[:, 0])
unique_ys = np.unique(xys[:, 1])
num_xs = len(unique_xs)
num_ys = len(unique_ys)

X, Y = np.meshgrid(unique_xs, unique_ys)

true_density = true_density.reshape(num_ys, num_xs)
initial_density = initial_density.reshape(num_ys, num_xs)
learned_density = learned_density.reshape(num_ys, num_xs)

# Define the mosaic layout
layout = """
    AAABBB
    CCDDEE
"""

# Create the figure with the specified layout
fig, axs = plt.subplot_mosaic(layout, figsize=(18, 12), constrained_layout=True)

# Plotting - tuning curves and estimation difference
axs['A'].contourf(X, Y, stcs.reshape(X.shape), levels=100, cmap='viridis')
axs['A'].scatter(pss[:, 0], pss[:, 1], color='red', label='Preferred Stimuli', s=2)
axs['A'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))
axs['A'].legend(loc='upper right')

axs['B'].contourf(X, Y, stcerrs.reshape(X.shape), levels=100, cmap='magma')
axs['B'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))

# Plotting - true, initial, and learned densities
for key, density, title in zip(['C', 'D', 'E'], [true_density, initial_density, learned_density], ['True Density', 'Initial Density', 'Learned Density']):
    cont = axs[key].contourf(X, Y, density, levels=20, cmap='viridis')
    # fig.colorbar(cont, ax=axs[key])
    axs[key].set_title(title)

plt.tight_layout()

plot_file_path = get_plot_path("population-code-2d-gaussian.png")

# Save the plot to a file
plt.savefig(plot_file_path)

