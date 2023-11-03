import matplotlib.pyplot as plt
import numpy as np
import json

# Local
from common import get_result_path, get_plot_path

# Replace this with the actual path to your JSON file
json_file_path = get_result_path('von-mises-mixture.json')

with open(json_file_path, 'r') as f:
    data = json.load(f)

# Extract data for plotting
samples = np.array(data['samples'])
true_confidence = np.array(data['true-confidence'])
em_confidence = np.array(data['em-confidence'])
ml_confidence = np.array(data['ml-confidence'])
sem_confidence = np.array(data['sem-confidence'])
sml_confidence = np.array(data['sml-confidence'])
cross_entropy = np.array(data['cross-entropy'])

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot confidence intervals
axes[0].scatter(samples[:, 0], samples[:, 1], c='black', label='Samples')
for i,crcl in enumerate(true_confidence):
    if i == 0:
        axes[0].plot(crcl[:, 0], crcl[:, 1], label='True', color='black')
    else:
        axes[0].plot(crcl[:, 0], crcl[:, 1], color='black')
for i,crcl in enumerate(em_confidence):
    if i == 0:
        axes[0].plot(crcl[:, 0], crcl[:, 1], label='EM')
    else:
        axes[0].plot(crcl[:, 0], crcl[:, 1])
for i,crcl in enumerate(ml_confidence):
    if i == 0:
        axes[0].plot(crcl[:, 0], crcl[:, 1], label='ML')
    else:
        axes[0].plot(crcl[:, 0], crcl[:, 1])
for i,crcl in enumerate(sem_confidence):
    if i == 0:
        axes[0].plot(crcl[:, 0], crcl[:, 1], label='SEM')
    else:
        axes[0].plot(crcl[:, 0], crcl[:, 1])
for i,crcl in enumerate(sml_confidence):
    if i == 0:
        axes[0].plot(crcl[:, 0], crcl[:, 1], label='SML')
    else:
        axes[0].plot(crcl[:, 0], crcl[:, 1])
axes[0].legend()
axes[0].set_title('Confidence Intervals')
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[0].set_xlim([0, 2*np.pi])
axes[0].set_ylim([0, 2*np.pi])
axes[0].set_xticks([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
axes[0].set_xticklabels(['0', '0.5π', 'π', '1.5π', '2π'])
axes[0].set_yticks([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
axes[0].set_yticklabels(['0', '0.5π', 'π', '1.5π', '2π'])

# Plot cross-entropy descent
axes[1].set_title('Cross-Entropy Descent')
axes[1].set_xlabel('Iterations')
axes[1].set_ylabel('Log-Likelihood')
epochs = np.arange(len(cross_entropy))
axes[1].plot(epochs, cross_entropy[:, 0], label='True', color='black', linewidth=2)
axes[1].plot(epochs, cross_entropy[:, 1], label='EM', linewidth=2)
axes[1].plot(epochs, cross_entropy[:, 2], label='ML', linewidth=2)
axes[1].plot(epochs, cross_entropy[:, 3], label='SEM', linewidth=2)
axes[1].plot(epochs, cross_entropy[:, 4], label='SML', linewidth=2)
axes[1].legend()
axes[1].set_xlim([0, max(epochs)])
axes[1].set_ylim([min(cross_entropy.flatten()), max(cross_entropy.flatten())])

# Adjust layout
plt.tight_layout()

plot_file_path = get_plot_path("von-mises-mixture.png")

# Save the plot to a file
plt.savefig(plot_file_path)
