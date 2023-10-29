import json
import matplotlib.pyplot as plt
import numpy as np

# Local
from common import get_result_path, get_plot_path

# File names
dirichlet_filename = "multivariate-dirichlet.json"
mvn_filename = "multivariate-normal.json"

# Load Dirichlet data
with open(get_result_path(dirichlet_filename), "r") as f:
    dirichlet_data = json.load(f)

# Load MVN data
with open(get_result_path(mvn_filename), "r") as f:
    mvn_data = json.load(f)

# Dirichlet data extraction
dxs = np.array(dirichlet_data["xrange"])
dys = np.array(dirichlet_data["yrange"])
dtrue_density = np.array(dirichlet_data["true-density"])
dlearned_density = np.array(dirichlet_data["learned-density"])
dlog_likelihood = dirichlet_data["log-likelihood"]
dsamples = np.array(dirichlet_data["samples"])

# MVN data extraction
mxs = np.array(mvn_data["xrange"])
mys = np.array(mvn_data["yrange"])
mtrue_density = np.array(mvn_data["true-density"])
msource_fit_density = np.array(mvn_data["source-fit-density"])
mnatural_fit_density = np.array(mvn_data["natural-fit-density"])
msamples = np.array(mvn_data["samples"])

plt.figure(figsize=(12, 6))

# Dirichlet plot
plt.subplot(1, 2, 1)
dirichlet_sample_scatter = plt.scatter(
    dsamples[:, 0], dsamples[:, 1], s=10, color="black", alpha=0.5
)

dirichlet_true_contour = plt.contour(
    dxs, dys, dtrue_density.reshape(len(dxs), len(dys)), colors="black"
)
dirichlet_learned_contour = plt.contour(
    dxs, dys, dlearned_density.reshape(len(dxs), len(dys)), colors="red"
)

# Extract legend elements from contour plots
h0 = dirichlet_sample_scatter
h1, _ = dirichlet_true_contour.legend_elements()
h2, _ = dirichlet_learned_contour.legend_elements()

# Legend
ax2 = plt.gca().inset_axes([0.55, 0.55, 0.4, 0.4])
ax2.plot(dlog_likelihood, color="black")
ax2.set_title("Dirichlet Log-Likelihood Training")
ax2.set_xlabel("Epochs")
ax2.set_ylabel("Log-Likelihood")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Dirichlet Distributions")
# legend bottom left
plt.legend(loc=3)
plt.legend(
    [h0, h1[0], h2[0]],
    ["Samples", "Dirichlet True Density", "Dirichlet Learned Density"],
    loc=3,
)

# MVN plot
plt.subplot(1, 2, 2)
mvn_sample_scatter = plt.scatter(
    msamples[:, 0], msamples[:, 1], s=10, color="black", label="MVN Samples"
)
mvn_true_contour = plt.contour(
    mxs,
    mys,
    mtrue_density.reshape(len(mxs), len(mys)),
    colors="black",
)
mvn_source_contour = plt.contour(
    mxs,
    mys,
    msource_fit_density.reshape(len(mxs), len(mys)),
    colors="blue",
)
mvn_natural_contour = plt.contour(
    mxs,
    mys,
    mnatural_fit_density.reshape(len(mxs), len(mys)),
    colors="red",
    linestyles="dashed",
)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Multivariate Normal Distributions")

g0 = mvn_sample_scatter
g1, _ = mvn_true_contour.legend_elements()
g2, _ = mvn_source_contour.legend_elements()
g3, _ = mvn_natural_contour.legend_elements()
plt.legend(
    [g0, g1[0], g2[0], g3[0]],
    [
        "Samples",
        "MVN True Density",
        "MVN Source Fit Density",
        "MVN Natural Fit Density",
    ],
)

plt.tight_layout()

plot_file_path = get_plot_path("multivariate-distributions.png")

# Save the plot to a file
plt.savefig(plot_file_path)
