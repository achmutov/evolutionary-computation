# /// script
# dependencies = [
#   "pandas",
#   "matplotlib",
#   "seaborn",
# ]
# ///

import sys

if len(sys.argv) == 1:
    cmd = sys.argv[0]
    print(f"usage: {cmd} INPUT_CSV")
    sys.exit(1)
inp_path = sys.argv[1]

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

sns.set_theme()

inp = pd.read_csv(inp_path)

for (instance, method), row in inp.groupby(["instance", "method"]):
    instance_data = pd.read_csv(instance, delimiter=";", header=None)
    instance_data.columns = ["x", "y", "cost"]

    indices_str = row["best_solution"].values[0]
    indices = [int(i) for i in indices_str.split(" ")]

    sns.lineplot(instance_data.iloc[indices], x='x', y='y', c='lightgrey', sort=False, zorder=1)
    sns.lineplot(instance_data.iloc[[indices[0], indices[-1]]], x='x', y='y', c='lightgrey', sort=False, zorder=1)
    sns.scatterplot(instance_data, x='x', y='y', hue="cost", zorder=2)
    plt.show()
