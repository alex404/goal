import json
import matplotlib.pyplot as plt

# Local
from common import get_result_path, get_plot_path

# Load data from the JSON file
with open(get_result_path('factor-analysis.json'), 'r') as file:
    data = json.load(file)

# Extract the data for plotting
conjugation_curve_data = data['conjugation-curve']
kalman_filter_data = data['kalman-filter']

# Prepare the data for plotting
# Unzip conjugation curve data
xsmps, ys, yhts = zip(*conjugation_curve_data)

# Transpose Kalman filter data
xpth0, zpth0, mus, sds, mus1, sds1 = zip(*kalman_filter_data)

# Set up the plot with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot the conjugation curve on the first subplot (ax1)
ax1.plot(xsmps, ys, 'k-', lw=4, label='Conjugation Curve')
ax1.plot(xsmps, yhts, 'r--', lw=4, label='Fit')
ax1.set_xlabel('x')
ax1.set_ylabel('Potential')
ax1.set_title('Conjugation Curve')
ax1.legend()

# Plot the Kalman filter results on the second subplot (ax2)
time = range(len(xpth0))
ax2.plot(time, xpth0, 'k-', lw=4, label='Latent State')
ax2.scatter(time, zpth0, color='black', label='Observations')
ax2.plot(time, mus, 'r-', lw=4, label='Optimal')
ax2.fill_between(time, [m + s for m, s in zip(mus, sds)], [m - s for m, s in zip(mus, sds)], color='red', alpha=0.5)
ax2.plot(time, mus1, 'b-', lw=4, label='Learned')
ax2.fill_between(time, [m + s for m, s in zip(mus1, sds1)], [m - s for m, s in zip(mus1, sds1)], color='blue', alpha=0.5)
ax2.set_xlabel('Time')
ax2.set_ylabel('Latent State')
ax2.set_title('Kalman Filter')
ax2.legend()

# Adjust the layout
plt.tight_layout()

plot_file_path = get_plot_path("kalman-filter.png")

# Save the plot to a file
plt.savefig(plot_file_path)


