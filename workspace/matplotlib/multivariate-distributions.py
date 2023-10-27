import json
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def plot_combined(directory):

    dirichlet_filename = "multivariate-dirichlet.json"
    mvn_filename = "multivariate-normal.json"


    # Load Dirichlet data
    with open(os.path.join(directory, dirichlet_filename), "r") as f:
        dirichlet_data = json.load(f)

    # Load MVN data
    with open(os.path.join(directory, mvn_filename), "r") as f:
        mvn_data = json.load(f)

    # Dirichlet data extraction
    dxs = np.array(dirichlet_data['xrange'])
    dys = np.array(dirichlet_data['yrange'])
    dtrue_density = np.array(dirichlet_data['true-density'])
    dlearned_density = np.array(dirichlet_data['learned-density'])
    dlog_likelihood = dirichlet_data['log-likelihood']
    dsamples = np.array(dirichlet_data['samples'])

    # MVN data extraction
    mxs = np.array(mvn_data['xrange'])
    mys = np.array(mvn_data['yrange'])
    mtrue_density = np.array(mvn_data['true-density'])
    msource_fit_density = np.array(mvn_data['source-fit-density'])
    mnatural_fit_density = np.array(mvn_data['natural-fit-density'])
    msamples = np.array(mvn_data['samples'])

    plt.figure(figsize=(15, 7))

    # Dirichlet plot
    plt.subplot(1, 2, 1)
    plt.scatter(dsamples[:, 0], dsamples[:, 1], s=10, color='red', label='Dirichlet Samples')
    plt.contour(dxs, dys, dtrue_density.reshape(len(dxs), len(dys)), colors='blue', label='Dirichlet True Density')
    plt.contour(dxs, dys, dlearned_density.reshape(len(dxs), len(dys)), colors='green', label='Dirichlet Learned Density')
    ax2 = plt.gca().inset_axes([0.5, 0.5, 0.4, 0.4])
    ax2.plot(dlog_likelihood)
    ax2.set_title('Dirichlet Log-Likelihood Training')
    ax2.set_xlabel('Epochs')
    ax2.set_ylabel('Log-Likelihood')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Dirichlet Distributions')
    plt.legend()

    # MVN plot
    plt.subplot(1, 2, 2)
    plt.scatter(msamples[:, 0], msamples[:, 1], s=10, color='red', label='MVN Samples')
    plt.contour(mxs, mys, mtrue_density.reshape(len(mxs), len(mys)), colors='blue', label='MVN True Density')
    plt.contour(mxs, mys, msource_fit_density.reshape(len(mxs), len(mys)), colors='green', label='MVN Source Fit Density')
    plt.contour(mxs, mys, mnatural_fit_density.reshape(len(mxs), len(mys)), colors='yellow', label='MVN Natural Fit Density')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Multivariate Normal Distributions')
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    plot_combined(sys.argv[1])

