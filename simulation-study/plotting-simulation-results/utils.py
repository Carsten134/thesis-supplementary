import matplotlib.pyplot as plt
import itertools
import pandas as pd

def plot_power_result(df):
    """
    Plots power results from a DataFrame.
    
    Parameters:
    df (pd.DataFrame): Must contain columns:
                       - 'decision': The rejection rate (y-axis)
                       - 'delta': The parameter delta (x-axis)
                       - 'n': The grid size identifier
    """
    
    # --- 1. Global Configuration (Academic Standard) ---
    plt.rcParams.update({
        "text.usetex": True,                # Use LaTeX rendering
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "legend.fontsize": 10,
        "lines.linewidth": 1.5,
        "lines.markersize": 6
    })

    # --- 2. Setup Figure ---
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=300)
    
    # Cyclers for visual distinction (B&W safe)
    markers = itertools.cycle(['o', 's', '^', 'D', 'v', 'x'])
    linestyles = itertools.cycle(['-', '--', '-.', ':', (0, (3, 1, 1, 1))])
    colors = itertools.cycle(["#333333", "#666666", "#999999", "#8B8B8B"])

    # --- 3. Grouping and Plotting ---
    # Get unique N values and sort them so the legend is ordered (20, 30, 40...)
    unique_ns = sorted(df["n"].unique())

    for n in unique_ns:
        # Filter and sort by delta to ensure lines are drawn correctly
        subset = df[df["n"] == n].sort_values("delta")
        
        # Format label: 20 -> 20 x 20
        label_str = r"${0} \times {0}$".format(int(n))
        
        ax.plot(
            subset["delta"], 
            subset["decision"], 
            label=label_str,
            marker=next(markers),
            linestyle=next(linestyles),
            color=next(colors),
            alpha=0.9
        )

    # --- 4. Labels and Titles ---
    # ax.set_title(r"$N = 200$ (Power Function)", pad=15)
    ax.set_ylabel("Approximated Power")
    ax.set_xlabel(r"$\delta$")

    # --- 5. Visual Polish ---
    ax.grid(True, which='major', linestyle='-', alpha=0.3, color='gray')
    ax.legend(title="Sample Size", frameon=True, fancybox=False, edgecolor='black', framealpha=1)
    
    # Optional: Set limits if data is strictly 0-1
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.show()