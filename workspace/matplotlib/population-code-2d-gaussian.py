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
tc_ys = np.array(data["tuning-curve-ys"])
y_samples = np.array(data["y-samples"])
stcs = np.array(data["sum-of-tuning-curves"])
stcerrs = np.array(data["estimation-difference"])
pss = np.array(data["preferred-stimuli"])
density_ys = np.array(data["density-ys"])
# true_density = np.array(data["true-density"])
initial_density = np.array(data["initial-density"])
learned_density = np.array(data["learned-density"])
regmn, regmx = data["regression-bounds"]

# Prepare the grid
unique_tc_y1s = np.unique(tc_ys[:, 0])
unique_tc_y2s = np.unique(tc_ys[:, 1])
unique_density_y1s = np.unique(density_ys[:, 0])
unique_density_y2s = np.unique(density_ys[:, 1])
num_y1s = len(unique_tc_y1s)
num_y2s = len(unique_tc_y2s)

tc_Y1, tc_Y2 = np.meshgrid(unique_tc_y1s, unique_tc_y2s)
density_Y1, density_Y2 = np.meshgrid(unique_density_y1s, unique_density_y2s)

# true_density = true_density.reshape(num_y2s, num_y1s)
initial_density = initial_density.reshape(num_y2s, num_y1s)
learned_density = learned_density.reshape(num_y2s, num_y1s)

spike_counts = np.array(data["x-samples"][0])
# Define the mosaic layout
layout = """
    AAABBB
    CCDDEE
"""

# Create the figure with the specified layout
fig, axs = plt.subplot_mosaic(layout, figsize=(18, 12), constrained_layout=True)

# Plotting - tuning curves and estimation difference
axs['A'].contourf(tc_Y1, tc_Y2, stcs.reshape(tc_Y1.shape), levels=100, cmap='viridis')
axs['A'].scatter(pss[:, 0], pss[:, 1], color='red', label='Preferred Stimuli', s=2)
axs['A'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))
axs['A'].legend(loc='upper right')

axs['B'].contourf(tc_Y1, tc_Y2, stcerrs.reshape(tc_Y1.shape), levels=100, cmap='magma')
axs['B'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))

# Scatter plot of the samples
axs['C'].scatter(y_samples[:, 0], y_samples[:, 1], s=2)
axs['C'].scatter(y_samples[0, 0], y_samples[0, 1], s=40, edgecolor='black', facecolor='red', label='Sample Point')
# Range based on the tuning curve bounds
axs['C'].add_patch(Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none'))

sizes = 15 * (1 + 6*spike_counts / np.max(spike_counts))
axs['D'].scatter(pss[:, 0], pss[:, 1], c=spike_counts, s=sizes, cmap='viridis')
# Add colorbar
fig.colorbar(axs['D'].collections[0], ax=axs['D'], orientation='horizontal')


# Plotting - true, initial, and learned densities
cont = axs['E'].contourf(density_Y1, density_Y2, learned_density, levels=20, cmap='viridis')
# fig.colorbar(cont, ax=axs[key])
axs['E'].set_title("Learned Prior Density")

plot_file_path = get_plot_path("population-code-2d-gaussian.png")

# Save the plot to a file
plt.savefig(plot_file_path)

