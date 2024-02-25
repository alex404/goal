### Imports ###

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Local
from common import get_result_path, get_plot_path

### Helper Functions ###

# Simplex plot helper
def plot_simplex(ax, weights, colors, **kwargs):
    # Triangle vertices for a 2D simplex representing a 3-component mixture
    vertices = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    triangle = Polygon(vertices, edgecolor='k', fill=False)
    ax.add_patch(triangle)

    # Plot the vertices with the corresponding colors
    for vertex, color in zip(vertices, colors):
        ax.plot(vertex[0], vertex[1], 'o', color=color, ms=10)

    # Plot the current weights as a point in the simplex
    barycentric = np.dot(vertices.T, weights)
    ax.plot(barycentric[0], barycentric[1], 'o', ms=10, color=weights)

    # Configure plot
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.axis('off')

### Initialization ###

# Load the JSON data
with open(get_result_path('conjugated-harmoniums.json'), 'r') as file:
    data = json.load(file)


# Layout for the plots: three rows with two columns
layout = """
AB
CD
EF
"""

# Create the plots using mosaic
fig, axs = plt.subplot_mosaic(layout, figsize=(14, 12), constrained_layout=True)

### Factor Analysis ###

# FA model parameters
W = np.array([2, 1])
mu0 = np.array([1, -1])
sigmas = np.array([np.sqrt(3), np.sqrt(2)])  # Standard deviations

# The mean of the FA model (for z=0)
mean = W * 0 + mu0

# Plot for FA observable densities (contour plot)
observable_x, observable_y, observable_density = zip(*data["fa-observable-densities"])
xi, yi = np.linspace(min(observable_x), max(observable_x), 100), np.linspace(min(observable_y), max(observable_y), 100)
zi = np.reshape(observable_density, (100, 100))
cs = axs['A'].contour(xi, yi, zi, cmap='viridis')
plt.clabel(cs, inline=True)
axs['A'].set_title('Factor Analysis Observable Density')

# cbar_ax = fig.add_axes([0.9, 0.55, 0.01, 0.3])  # x, y, width, height in figure coordinate
# Plot for FA latent densities
latent_x, latent_density = zip(*data["fa-latent-densities"])
axs['B'].plot(latent_x, latent_density, color='k')
axs['B'].set_title('Factor Analysis Prior')

# Set x and y axis ticks to show the mean and +/- 1 standard deviation
# Note: We'll use symbolic representation here
axs['A'].set_xticks([mean[0] - sigmas[0], mean[0], mean[0] + sigmas[0]])
axs['A'].set_yticks([mean[1] - sigmas[1], mean[1], mean[1] + sigmas[1]])
axs['A'].set_ylim([mean[1] - 1.75*sigmas[1], mean[1] + 1.75*sigmas[1]])

# Symbolic labels
axs['A'].set_xticklabels(
    [r'$\mu_x + w_x z - \sigma_x$' , r'$\mu_x + w_x z$', r'$\mu_x + w_x z + \sigma_x$'])
axs['A'].set_yticklabels(
    [r'$\mu_y + w_y z - \sigma_y$', r'$\mu_y + w_y z$', r'$\mu_y + w_y z + \sigma_y$']
    , va='center', rotation=90)### Mixture Model ###

# Define colors for the mixture components
colors = ['red', 'green', 'blue']  # Colors for the mixture components

# Extract mixture weights from the data
mixture_weights = np.array(data["mixture-weights"])


# Plot for mixture model densities (panel 'C')
mix_samples = data["mixture-samples"]
mix_densities = data["mixture-densities"]
axs['C'].plot(mix_samples, mix_densities, label='Mixture Density', color='black')

# Plot the individual component densities
for comp_density, color in zip([data["mixture-component-1-densities"], 
                                data["mixture-component-2-densities"], 
                                data["mixture-component-3-densities"]], 
                               colors):
    axs['C'].plot(mix_samples, comp_density, label=f'Component Density', color=color)
axs['C'].legend()

# Plot the simplex with the mixture weights
plot_simplex(axs['D'], mixture_weights, colors)
axs['C'].set_title('Mixture and Component Densities')
axs['D'].set_title('Mixture Weights')

### Population Coding ###

# Assuming axs['F'] is the axis for the Von Mises density plot
axs['E'].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
axs['E'].set_xticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
axs['E'].set_xlabel('Angle (radians)')

# Plot for Neuron Responses (tuning curves)
for neuron_response in data["neuron-responses"]:
    axs['E'].plot(data["stimulus-samples"], neuron_response, color='grey', alpha=0.5)

# Assuming axs['F'] is the axis for the Von Mises density plot
axs['F'].set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
axs['F'].set_xticklabels(['0', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
axs['F'].set_xlabel('Angle (radians)')

# Plot for Von Mises Density (prior)
axs['F'].plot(data["stimulus-samples"], data["von-mises-density"], color='k')

# Customize axis labels for the tuning curves and Von Mises prior
axs['E'].set_xlabel('Stimulus')
axs['E'].set_ylabel('Response')
axs['F'].set_xlabel('Stimulus')
axs['E'].set_title('Tuning Curves')
axs['F'].set_title('Von Mises Prior')

# Panel labels
# for label, ax in axs.items():
#     ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top')

# Save and show the plot
plt.savefig(get_plot_path("conjugated-harmoniums.png"))

