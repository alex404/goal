import os
import matplotlib.pyplot as plt

# Determine the directory containing the currently executing script
_script_dir = os.path.dirname(os.path.realpath(__file__))

# Construct paths to the 'results' and 'plots' directories
rsltsdr = os.path.join(_script_dir, "..", "results")
pltsdr = os.path.join(_script_dir, "..", "plots")

plt.style.use(os.path.join(_script_dir, "default.mplstyle"))

# Ensure the plots directory exists
if not os.path.exists(pltsdr):
    os.makedirs(pltsdr)


def get_result_path(filename):
    """Returns the full path for a given filename in the results directory."""
    return os.path.join(rsltsdr, filename)


def get_plot_path(filename):
    """Returns the full path for a given filename in the plots directory."""
    return os.path.join(pltsdr, filename)
