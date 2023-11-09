import json
import matplotlib.pyplot as plt

# Local
from common import get_result_path, get_plot_path

# Load data from the JSON file
with open(get_result_path('factor-analysis.json'), 'r') as file:
    data = json.load(file)

# Extract the data for plotting
standard_log_likelihood = data['standard-log-likelihood']
ef_log_likelihood = data['ef-log-likelihood']
data_correlations = data['data-correlations']
standard_factor_correlations = data['standard-factor-analysis-correlations']
ef_factor_correlations = data['ef-factor-analysis-correlations']

# Set up the plot with two subplots: one for learning trajectory, one for correlations
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot the learning trajectory on the left subplot (ax1)
ax1.plot(range(len(standard_log_likelihood)), standard_log_likelihood, label='Standard FA Log Likelihood')
ax1.plot(range(len(ef_log_likelihood)), ef_log_likelihood, label='EF FA Log Likelihood')
ax1.set_title('Learning Trajectory')
ax1.set_xlabel('Iteration')
ax1.set_ylabel('Log Likelihood')
ax1.legend()

# Plot the correlations on the right subplot (ax2)
ax2.scatter(data_correlations, standard_factor_correlations, c='red', label='Standard FA Correlations')
ax2.scatter(data_correlations, ef_factor_correlations, c='blue', label='EF FA Correlations')
ax2.plot([-1, 1], [-1, 1], 'k--', label='y = x')  # Line y = x for reference
ax2.set_title('Correlations Scatter')
ax2.set_xlabel('Data Correlations')
ax2.set_ylabel('FA Correlations')
ax2.legend()
ax2.axis('equal')  # Set the same scale for x and y axis

# Adjust the layout
plt.tight_layout()

plot_file_path = get_plot_path("factor-analysis.png")

# Save the plot to a file
plt.savefig(plot_file_path)

