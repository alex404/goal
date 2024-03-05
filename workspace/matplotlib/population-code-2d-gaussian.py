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

xys = np.array(data["xys"])  # Assuming this is a list of [x, y] coordinates
stcs = np.array(data["sum-of-tuning-curves"])  # A list of lists, each sublist is a tuning curve across xys
stcerrs = np.array(data["estimation-difference"])  # A list of lists, each sublist is a tuning curve across xys
regression_bounds = data["regression-bounds"]
pss = data["preferred-stimuli"]
regmn, regmx = regression_bounds

# Reshape xys to facilitate contour plotting
x_vals = sorted(set([xy[0] for xy in xys]))
y_vals = sorted(set([xy[1] for xy in xys]))

X, Y = np.meshgrid(x_vals, y_vals)
# Create a figure for plotting
fig, axes = plt.subplots(nrows=2, figsize=(8, 12))

### Set global configurations for both plots
lvls = 100

### Plot the sum of tuning curves
contourtc = axes[0].contourf(X, Y, stcs.reshape(X.shape), levels=lvls, cmap='viridis')
reg_box1 = Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none')
sct = axes[0].scatter([ps[0] for ps in pss], [ps[1] for ps in pss], color='red', label='Preferred Stimuli', s=2)
axes[0].add_patch(reg_box1)

cbartc = fig.colorbar(contourtc, ax=axes[0])
cbartc.set_label('Firing Rate')
axes[0].set_xlabel('X Coordinate')
axes[0].set_ylabel('Y Coordinate')
axes[0].set_title('Sum of 2D Gaussian Tuning Curves')

scatter_legend = Line2D([0], [0], linestyle="none", color='red', marker='o') 
reg_legend = Line2D([0], [0], color='black', lw=2)
axes[0].legend([scatter_legend, reg_legend], ['Preferred Stimuli', 'Regression Bounds'], loc='upper right')

### Plot the estimation difference
contourerr = axes[1].contourf(X, Y, stcerrs.reshape(X.shape), levels=lvls, cmap='magma')
reg_box2 = Rectangle((regmn, regmn), regmx-regmn, regmx-regmn, linewidth=2, edgecolor='black', facecolor='none')
axes[1].add_patch(reg_box2)

cbarerr = fig.colorbar(contourerr, ax=axes[1])
cbarerr.set_label('Estimation Residual')
axes[1].set_xlabel('X Coordinate')
axes[1].set_ylabel('Y Coordinate')
axes[1].set_title('Estimation Residual')

# Adjust the layout
plt.tight_layout()

plot_file_path = get_plot_path("population-code-2d-gaussian.png")

# Save the plot to a file
plt.savefig(plot_file_path)

