import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.patches import Rectangle

# Local
from common import get_result_path, get_plot_path

### Load and prepare the data
with open(get_result_path("population-code-2d-gaussian.json"), "r") as file:
    data = json.load(file)

# Extracting the necessary elements from the data
tc_xys = np.array(data["tuning-curve-xys"])
stcs = np.array(data["sum-of-tuning-curves"])
stcerrs = np.array(data["estimation-difference"])
pss = np.array(data["preferred-stimuli"])
density_xys = np.array(data["density-xys"])
true_density = np.array(data["true-density"])
initial_density = np.array(data["initial-density"])
learned_density = np.array(data["learned-density"])
regmn, regmx = data["regression-bounds"]

# Prepare the grid
unique_tc_xs = np.unique(tc_xys[:, 0])
unique_tc_ys = np.unique(tc_xys[:, 1])
unique_density_xs = np.unique(density_xys[:, 0])
unique_density_ys = np.unique(density_xys[:, 1])
num_xs = len(unique_tc_xs)
num_ys = len(unique_tc_ys)

tc_X, tc_Y = np.meshgrid(unique_tc_xs, unique_tc_ys)
density_X, density_Y = np.meshgrid(unique_density_xs, unique_density_ys)

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
axs['A'].contourf(tc_X, tc_Y, stcs.reshape(tc_X.shape), levels=100, cmap='viridis')
axs['A'].scatter(pss[:, 0], pss[:, 1], color='red', label='Preferred Stimuli', s=2)
axs['A'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))
axs['A'].legend(loc='upper right')

axs['B'].contourf(tc_X, tc_Y, stcerrs.reshape(tc_X.shape), levels=100, cmap='magma')
axs['B'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))

# Plotting - true, initial, and learned densities
for key, density, title in zip(['C', 'D', 'E'], [true_density, initial_density, learned_density], ['True Prior Density', 'Initial Prior Density', 'Learned Prior Density']):
    cont = axs[key].contourf(density_X, density_Y, density, levels=20, cmap='viridis')
    # fig.colorbar(cont, ax=axs[key])
    axs[key].set_title(title)

plot_file_path = get_plot_path("population-code-2d-gaussian.png")

# Save the plot to a file
plt.savefig(plot_file_path)

