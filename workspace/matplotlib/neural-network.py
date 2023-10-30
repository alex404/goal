import matplotlib.pyplot as plt
import json

# Local
from common import get_result_path, get_plot_path

# Load the data
with open(get_result_path("neural-network.json"), "r") as f:
    data = json.load(f)

# Extract the x and y values of the samples
x_samples = [pair[0] for pair in data["xys"]]
y_samples = [pair[1] for pair in data["xys"]]

# Create a 2x1 subplot layout
fig, axs = plt.subplots(2, 1, figsize=(8, 8))

# Top panel: samples, true function, learned regressions
axs[0].scatter(x_samples, y_samples, label="Samples", color="black", marker='x')
axs[0].plot(data["inputs"], data["true-outputs"], label="True Function", color="gray", linewidth=2)
axs[0].plot(data["inputs"], data["gradient-descent-outputs"], label="Gradient Descent", color="blue", linewidth=2)
axs[0].plot(data["inputs"], data["momentum-outputs"], label="Momentum", color="green", linewidth=2)
axs[0].plot(data["inputs"], data["adam-outputs"], label="Adam", color="red", linewidth=2)
axs[0].legend()
axs[0].set_xlabel("Input")
axs[0].set_ylabel("Output")
axs[0].set_title("Samples, True Function, and Learned Regressions")

# Bottom panel: training histories
epochs = list(range(len(data["gradient-descent-training"])))
axs[1].plot(epochs, data["gradient-descent-training"], label="Gradient Descent", color="blue")
axs[1].plot(epochs, data["momentum-training"], label="Momentum", color="green")
axs[1].plot(epochs, data["adam-training"], label="Adam", color="red")
axs[1].legend()
axs[1].set_xlabel("Epochs")
axs[1].set_ylabel("Log-Likelihood")
axs[1].set_title("Training Histories")

plt.tight_layout()

plot_file_path = get_plot_path("neural-network.png")

# Save the plot to a file
plt.savefig(plot_file_path)


