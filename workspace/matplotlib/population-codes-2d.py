import matplotlib.pyplot as plt
import numpy as np
import json

# Local
from common import get_result_path, get_plot_path

# Load the data from the JSON file
with open(get_result_path("population-codes-2d.json"), "r") as file:
    data = json.load(file)

xys = np.array(data["xys"])  # Assuming this is a list of [x, y] coordinates
stcs = np.array(data["sum-of-tuning-curves"])  # A list of lists, each sublist is a tuning curve across xys

# Reshape xys to facilitate contour plotting
x_vals = np.linspace(-3, 3, 200)  # Adjust according to the actual range and number of points used
y_vals = np.linspace(-3, 3, 200)  # Adjust accordingly
X, Y = np.meshgrid(x_vals, y_vals)

# Create a figure for plotting
fig, ax = plt.subplots(figsize=(8, 6))

# Plot contours for each tuning curve
# Note: Adjust the levels parameter as necessary to get a desired contour detail
Z = stcs.reshape(X.shape)  # Reshape the tuning curve data to match the grid
contour = ax.contourf(X, Y, Z, levels=20, cmap='viridis')  # Use viridis colormap, adjust as needed

# Add a color bar to indicate firing rate density
cbar = fig.colorbar(contour, ax=ax)
cbar.set_label('Firing Rate')

# Set labels and title
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_title('2D Gaussian Tuning Curves')

# Adjust the layout
plt.tight_layout()

plot_file_path = get_plot_path("population-codes-2d.png")

# Save the plot to a file
plt.savefig(plot_file_path)

