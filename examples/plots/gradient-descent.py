import json
import matplotlib.pyplot as plt
import numpy as np
import sys

def compute_step_sizes(points):
    """Compute step sizes based on point differences."""
    diff = np.diff(points, axis=0)
    return np.sqrt((diff ** 2).sum(axis=1))

def plot_steps(ax, points, color, label):
    """Plot the steps with size based on step magnitude."""
    step_sizes = compute_step_sizes(points)
    sizes = (step_sizes / step_sizes.max() * 100) + 20  # normalize and amplify sizes for visibility
    ax.scatter(points[:-1, 0], points[:-1, 1], c=color, s=2, label=label)
    ax.plot(points[:, 0], points[:, 1], color=color, alpha=0.5, lw=0.5)

# Load the data from the JSON file
json_pth = sys.argv[1]

with open(json_pth, 'r') as file:
    data = json.load(file)

# Extract the isosamples and reshape it for contour plotting
isosamples = np.array(data["isosamples"])
X = isosamples[:, 0].reshape(100, 100)
Y = isosamples[:, 1].reshape(100, 100)
Z = isosamples[:, 2].reshape(100, 100)

# Create contour plot
fig, ax = plt.subplots(figsize=(10, 7))
contour = ax.contourf(X, Y, Z, 20, cmap='viridis')
plt.colorbar(contour, ax=ax)

# Plot the Adam, Momentum, and Classic points
adam_points = np.array(data["adam"])
momentum_points = np.array(data["momentum"])
classic_points = np.array(data["gradient-descent"])

plot_steps(ax, classic_points, 'cyan', 'Classic')
plot_steps(ax, momentum_points, 'magenta', 'Momentum')
plot_steps(ax, adam_points, 'yellow', 'Adam')

ax.legend(loc='upper right')
ax.set_title("Gradient Descent Comparison")
plt.show()

