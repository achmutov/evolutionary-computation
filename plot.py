# /// script
# dependencies = [
#   "pandas",
#   "matplotlib",
#   "seaborn",
# ]
# ///

import sys
from pathlib import Path

if len(sys.argv) == 1:
    cmd = sys.argv[0]
    print(f"usage: {cmd} INPUT_CSV")
    sys.exit(1)
inp_path = sys.argv[1]

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

sns.set_theme()

# Read CSV - read all columns as strings first to avoid type inference issues
# Then convert numeric columns back to numeric types
inp = pd.read_csv(inp_path, dtype=str, keep_default_na=False)

# Convert numeric columns back to appropriate types (but keep best_solution as string)
numeric_cols = ['min', 'max', 'mean', 'duration_min', 'duration_max', 'duration_mean', 
                'ls_runs_min', 'ls_runs_max', 'ls_runs_mean']
for col in numeric_cols:
    if col in inp.columns:
        inp[col] = pd.to_numeric(inp[col], errors='coerce')

for (instance, method), row in inp.groupby(["instance", "method"]):
    instance_data = pd.read_csv(instance, delimiter=";", header=None)
    instance_data.columns = ["x", "y", "cost"]

    indices_str = str(row["best_solution"].values[0]).strip()
    if not indices_str or indices_str.lower() == "nan":
        print(f"Warning: No best_solution found for {instance} - {method}, skipping...")
        continue
    indices = [int(i) for i in indices_str.split(" ") if i.strip()]

    # Create complete cycle by adding first node at the end
    cycle_indices = indices + [indices[0]]
    
    # Extract coordinates as simple arrays
    x_coords = instance_data['x'].values
    y_coords = instance_data['y'].values
    costs = instance_data['cost'].values
    
    # Create path coordinates in correct order
    path_x = [x_coords[idx] for idx in cycle_indices]
    path_y = [y_coords[idx] for idx in cycle_indices]
    
    # Plot using matplotlib directly
    plt.plot(path_x, path_y, 'k-', linewidth=1, zorder=1)
    scatter = plt.scatter(x_coords, y_coords, c=costs, s=20, zorder=2, cmap='RdPu')
    
    # Add colorbar legend
    cbar = plt.colorbar(scatter)
    cbar.set_label('Node Cost', rotation=270, labelpad=15)
    
    # Add axis labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'{Path(instance).stem} - {method.replace("_", " ").title()} Solution')

    instance_path = Path(instance).with_suffix('')
    instance_plot_path = instance_path.with_stem(f"{instance_path.stem}-{method}").with_suffix(".png")

    plt.savefig(instance_plot_path, dpi=300)
    plt.clf()
