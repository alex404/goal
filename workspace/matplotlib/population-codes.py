import matplotlib.pyplot as plt
import numpy as np
import json

# Local
from common import get_result_path, get_plot_path

# Load the data from the JSON file
with open(get_result_path("population-codes.json"), "r") as file:
    data = json.load(file)

# Create a figure with two subplots
fig, axes = plt.subplots(nrows=2, figsize=(8, 8))

axes2 = []

# Set global configurations for both plots
for i, ax in enumerate(axes):
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks([0, np.pi/2, np.pi, 1.5*np.pi, 2*np.pi])
    ax.set_xticklabels(["0", "0.5π", "π", "1.5π", "2π"])
    ax.set_xlabel('Orientation')
    ax.set_ylabel('Rate')
    ax.set_ylim(0, 8)
    axes2.append(ax.twinx())
    axes2[i].set_ylabel('Total Rate')
    axes2[i].set_ylim(0, 30)
    axes2[i].spines['right'].set_visible(True)

# Plot for linear population
l1, = axes2[0].plot(data["xs"], data["linear-sum-of-tuning-curves"], lw=4, label="Sum of Tuning Curves")
l2, = axes2[0].plot(data["xs"], data["linear-regression"], lw=4, linestyle='--', label="Sinusoid Fit")
l3, = axes[0].plot(data["xs"], data["linear-tuning-curves"][0], lw=4, color='black', label="Tuning Curves")
for curve in data["linear-tuning-curves"][1:]:
    axes[0].plot(data["xs"], curve, lw=4, color='black')

# Plot for affine population
axes2[1].plot(data["xs"], data["affine-sum-of-tuning-curves"], lw=4)
axes2[1].plot(data["xs"], data["affine-regression"], lw=4, linestyle='--')
for curve in data["affine-tuning-curves"]:
    axes[1].plot(data["xs"], curve, lw=4, color='black')

# Unified legend
fig.legend(handles=[l1, l2, l3], loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.05))

# Adjust the layout
plt.tight_layout()

plot_file_path = get_plot_path("population-codes.png")

# Save the plot to a file
plt.savefig(plot_file_path)
