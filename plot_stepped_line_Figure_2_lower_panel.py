"""
Description:
Plots Figure 2's lower panel: the per–0.1% HbA1c-bin diabetes misclassification (multi separate lines for +1, -1), no histogram background from the ...misclasss_analysis.csv file

Author: Hugh Claridge
Start date: January 2025
Last edited date: 28th October 2026
Version: v1.0
Version history: Public version of the plotting code for manuscript Figure 2's lower panel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from cmcrameri import cm as crm

## Load in the analysis file for checking and track the filename for plot tracing
## By default this will run on the data file used in the paper
## If you have created a new data file then you will need to replace the below with the new file name
analysis_bins_file = "20260104_233659_n10000_s1000_4.5-9.4_misclass_analysis.csv" ## Data file for plotting

## Load per-bin summary
bin_stats = pd.read_csv(analysis_bins_file)

# Make the interval midpoints from the strings in the file
# Parse the left edge and add 0.05 to get the middle of it since bins are 0.1 wide
left_edges = bin_stats["HbA1c Bin"].str.extract(r"\[([\d\.]+),")[0].astype(float)
bin_midpoints = (left_edges + 0.05).tolist()

## Plotting
fig, ax1 = plt.subplots(figsize=(12, 5))                                  ## Create figure 12w x 5h inches, and x and y axes

## +1 category line
ax1.step(
    bin_midpoints, bin_stats["Mean1Above"],
    where="mid", linewidth=2.0, color=crm.lajolla(0.50),
    label="+1 category"
)
## +1 category SD ribbon
lower1above = np.clip(bin_stats["Mean1Above"] - bin_stats["SD1Above"], 0, 100)
upper1above = np.clip(bin_stats["Mean1Above"] + bin_stats["SD1Above"], 0, 100)
ax1.fill_between(
    bin_midpoints, lower1above, upper1above,
    alpha=0.32, step="mid", color=crm.lajolla(0.50) #, label="±1 SD (+1 category)"
)

## -1 category line
ax1.step(
    bin_midpoints, bin_stats["Mean1Below"],
    where="mid", linewidth=2.0, color=crm.lajolla(0.85),
    label="-1 category", zorder=6
)
## -1 category SD ribbon
lower1below = np.clip(bin_stats["Mean1Below"] - bin_stats["SD1Below"], 0, 100)
upper1below = np.clip(bin_stats["Mean1Below"] + bin_stats["SD1Below"], 0, 100)
ax1.fill_between(
    bin_midpoints, lower1below, upper1below,
    alpha=0.28, step="mid", color=crm.lajolla(0.85) #, label="±1 SD (-1 category)"
)

ax1.set_ylabel("Mean (%) of $t$-samples in a given distance bucket", color="black")
ax1.set_ylim(0, 100)                                                   ## Limit 0-100 as %
ax1.yaxis.set_major_locator(MultipleLocator(10))                       ## 0,10,20,...,100

## X-axis ticks/limits derived from the bin midpoints
xtick_min = np.floor(min(bin_midpoints) * 2) / 2
xtick_max = np.ceil(max(bin_midpoints) * 2) / 2
ax1.set_xticks(np.arange(xtick_min, xtick_max + 0.1, 0.5))   ## Add x axis markers every 0.5% HbA1c
ax1.set_xlim(xtick_min, xtick_max)                           ## Limit range
ax1.grid(True, linestyle='--', alpha=0.6)                    ## Faint dashed gridlines
ax1.set_xlabel("HbA1c% (0.1 increment bins)")

ax1.set_title("Distance from point estimate's diabetes category when converting HbA1c% to eAG")

# Combined legend
fig.legend(
    loc="center",                   ## Point in the legend where bbox anchors
    bbox_to_anchor=(0.835, 0.53),    ## x, y coordinates for legend
    ncol=1,
    frameon=True,
    borderaxespad=0.0
)

plt.annotate(f"Source: {analysis_bins_file}", xy=(0.01, 0.01), xycoords="figure fraction", fontsize=9, color="gray")  ## Put plot's data source filename in bottom left corner for traceability

# Save
out_plot = analysis_bins_file.replace("misclass_analysis.csv", "stepped_line_Figure_2_lower_panel.png")
plt.tight_layout()
plt.savefig(out_plot, dpi=400)
print(f"Plot saved as: {out_plot}")
plt.show()
