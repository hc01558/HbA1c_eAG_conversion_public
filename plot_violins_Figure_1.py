"""
Visualises three sources of uncertainty:
1. HbA1c rounding (0.1% increments)
2. HbA1c assay imprecision (jitter, as relative CV)
3. Conversion model uncertainty (t-distribution)

This version has one single point on the plot to show rounding uncertainty, rather than having the bars applied to every point - manually positioned point to enable it to be to scale.
Assay jitter is overlaid, it's not part of the simulation as otherwise couldn't disentangle the sources of uncertainty in the plot (as conversion uncertainty and assay jitter would be combined).
Diabetes bands are marked on the right hand y-axis with ticks rather than integrated as coloured shaded bands.

Glossary:
eAG: Estimated Average Glucose
HbA1c: Glycated haemoglobin/hemoglobin
NGSP: National Glycohemoglobin Standardization Program

Author: Hugh Claridge
Start date: January 2025
Last edited date: 28th October 2026
Version: v1.0
Version history: Public version of the plotting code for manuscript Figure 1
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm as crm
from core_formulae import ngsp_to_eag_formula
from core_converter_tools import round_half_up

filename = "20260104_233659_n10000_s1000_4.5-9.4.csv" ## Data file for use by the plotting code


df = pd.read_csv(filename)                                     ## Load in chosen csv file

# Only include rows with input HbA1c % between 4.5 and 9.5
df = df[(df["Original HbA1c"] >= 4.5) & (df["Original HbA1c"] <= 9.5)]

## Conversion model - extract all t-samples
t_cols = [c for c in df.columns if c.startswith("T-Sample-")]   ## Get the 1000 T_sample col names into a list
t_samples = df[t_cols]                                          ## Use the list to make a new df containing the values from those cols

## Make the plot
fig, ax = plt.subplots(figsize=(12, 5))                         ## Create figure 12w x 5h inches, and x and y axes

## Violin plots per 0.1% HbA1c bin
df["Hba1c_bin"] = df["Original HbA1c"].apply(lambda hba1c_raw: round_half_up(hba1c_raw, 1))         ## 0.1% bins based on Original HbA1c

bin_vals = []                                              ## List of arrays: all t-samples in each bin
bin_pos = []                                               ## x-positions (bin HbA1c values)
pi_lower = []                                              ## Per-bin median analytic 90% PI lower bound from spreadsheet
pi_upper = []                                              ## Per-bin median analytic 90% PI upper bound from spreadsheet

for hb in np.sort(df["Hba1c_bin"].unique()):                ## Loop over each 0.1% HbA1c bin value in ascending order
    rows = df["Hba1c_bin"] == hb                            ## Boolean mask to hide rows not in this bin (keep only rows where this is True)
    if not rows.any():                                      ## Safety skip, if no rows fall in this bin (shouldn’t happen)
        continue                                            ## Occasional yellow error on 'any' is PyCharm thinking it's a plain boolean when it isn't.
    arr = t_samples.loc[rows].values.ravel(order="C")       ## From those rows, take all T-Sample cols and flatten to 1D (row-major 'C' order) to create one long vector of all t-distribution samples for this bin.
    pi_lower.append(df.loc[rows, "Lower Bound"].median())   ## Per bin median of lower 90% PI limit calculated using Nathan's formula
    pi_upper.append(df.loc[rows, "Upper Bound"].median())   ## Per bin median of upper 90% PI limit calculated using Nathan's formula
    bin_vals.append(arr)                                    ## Store the 1D vector of samples for this bin - matplot's violin plotting needs this
    bin_pos.append(hb)                                      ## Store the bin’s HbA1c value x-position to place the violin

v = ax.violinplot(                                          ## Draw a set of violin plots
    bin_vals,                                               ## A list/sequence: per-bin 1-D arrays of samples (one array per HbA1c bin)
    positions=bin_pos,                                      ## X positions for each violin (the matching HbA1c bin centers)
    vert=True,                                              ## Draw vertically (x = bin, y = values)
    widths=0.08,                                            ## Violin width ~80% of a 0.1% HbA1c step
    showmeans=False, showmedians=True, showextrema=False    ## Show median line only
)                                                           ## Returns a container ‘v’ that holds the drawn violins for style tweaking

for b in v['bodies']:                                       ## Iterate each violin body (one PolyCollection per HbA1c bin)
    b.set_facecolor(crm.grayC(0.95))                        ## Fill colour of light grey - higher for ligher
    b.set_edgecolor(crm.grayC(0.35))                        ## Mid/dark grey edge
    b.set_linewidth(0.6)                                    ## Outer lines
    b.set_alpha(1.0)                                        ## Solid, reduce for transparency
    b.set_hatch("")                                         ## Texture switch out, looked confusing


## Add a legend handle for violins (Matplotlib doesn’t auto-legend them)
ax.plot([], [], lw=8, color=crm.grayC(0.85), alpha=0.25,
        label="Conversion uncertainty ($t$-samples, per 0.1 % bin)")

## Assay jitter here is overlaying error bars to the point estimates, no part of the simulation incorporates assay jitter so that the plot can present the sources of uncertainty separately
## Assay imprecision (jitter) +-2.3% standard uncertainty as half of the expanded uncertainty +-4.6% reported by Unal and Atikeler 2018, now going to be 2.3*k factor of 1.645 to get 90%
expanded_uncertainty = 0.038
## Assay jitter bounds per 0.1% bin
bin_hb = np.array(bin_pos)  ## HbA1c at bin centers
jitter_lower = ngsp_to_eag_formula(bin_hb * (1 - expanded_uncertainty))     ## Scale HbA1c down by 3.8%, convert to eAG
jitter_upper = ngsp_to_eag_formula(bin_hb * (1 + expanded_uncertainty))     ## Scale HbA1c up by 3.8%, convert to eAG

# Width calculation to use for both jitter lines and median tick lines
half_tick = 0.035  ## 35% of a 0.1% step gives total width (*2) of 0.07 for spacing between

for i, (x, ylo, yhi) in enumerate(zip(bin_pos, jitter_lower, jitter_upper), start=1):   ## Loop over the x position and two jitter bounds per bin
    ax.hlines([ylo, yhi], x - half_tick, x + half_tick,                     ## Draw lines at lower and upper jit bounds, left end cap to right end cap centred around bin value
              colors=crm.lajolla(0.6), linewidth=0.6, alpha=1.0,
              label=f"Assay imprecision (± {expanded_uncertainty*100:.1f} % expanded uncertainty, $k$ = 1.645)" if i == 1 else "_nolegend_")  ## Add label of first line only to legend, to remove need of full label on y axis


## P5–P95 around the median to show where middle 90% of the conversion samples lie
bin_medians = np.array([np.median(vals) for vals in bin_vals])         ## Violin median per bin
bin_p5 = np.array([np.percentile(vals, 5) for vals in bin_vals])    ## 5th percentile of the t-samples in the bin (5% of samples are below it)
bin_p95 = np.array([np.percentile(vals, 95) for vals in bin_vals])  ## 95th percentile of the t-samples in the bin (5% of samples are above it)

low_err = bin_medians - bin_p5
high_err = bin_p95 - bin_medians

ax.errorbar(
    bin_pos, bin_medians,
    yerr=[low_err, high_err],
    fmt="none",
    ecolor="black",
    elinewidth=1.2,
    capsize=3,
    label="P5–P95 ($t$-samples)"
)

# ## Overlay analytic 90% PI (Lower Bound/Upper Bound) - these are the ones calculated directly from the model
# pi_x_offset = 0.03  # small x offset so the PI sits a bit to the right of the P5–P95 bars
#
# for j, (x, lo, hi) in enumerate(zip(bin_pos, pi_lower, pi_upper)):
#     ax.vlines(
#         x + pi_x_offset, lo, hi,
#         color=crm.navia(0.6),
#         linewidth=0.8,
#         alpha=1.0,
#         zorder=3,
#         label="Analytic 90% PI" if j == 0 else "_nolegend_"
#     )

## Make the built-in violin median line the same x width as the jitter lines for clarity
median_segments = [np.array([[x - half_tick, m], [x + half_tick, m]])
                   for x, m in zip(bin_pos, bin_medians)]
v['cmedians'].set_segments(median_segments)                             ## Override default full-width segments to match the jitter lines
v['cmedians'].set_linewidth(0.6)                                        ## Median line y width (height)
v['cmedians'].set_color("black")
v['cmedians'].set_zorder(4)                                             ## Set behind the point-estimate dots (zorder=5)

## Point estimate reference point per bin by converting the bin HbA1c
point_estimate = ngsp_to_eag_formula(np.array(bin_pos))

ax.scatter(
    bin_pos, point_estimate,
    s=18, marker="D", facecolors="white", edgecolors=crm.grayC(0.20), linewidths=1.2, zorder=5,
    label="Point estimate"
)

## The following steps scale the X axis based on the inputs (for the varying ranges)
hba1c_min = df["Original HbA1c"].min()                                              ## Find lowest HbA1c input value in the dataset
hba1c_max = df["Original HbA1c"].max()                                              ## Find highest HbA1c input value in the dataset

## Round ticks to nearest 0.5 bounds
xtick_min = np.floor(hba1c_min * 2) / 2                           ## Round down to nearest 0.5 increment
xtick_max = np.ceil(hba1c_max * 2) / 2                            ## Round up to nearest 0.5 increment

ax.set_xticks(np.arange(xtick_min, xtick_max + 0.1, 0.5))               ## Add x axis markers every 0.5% HbA1c between bounds for visual help
ax.set_xlim(xtick_min, xtick_max)                                       ## Limit x axis range to the previously calculated min and max
ax.grid(axis="x", linestyle='--', alpha=0.6)                            ## Only x axis ticks so as to not confuse the two y axes

## Find the y axis limits
ymin, ymax = ax.get_ylim()
ax.set_ylim(0, ymax)  # Force bottom to 0

## Diabetes diagnostic thresholds in mg/dL on right hand axis with horizontal boundary lines - line labels may need changing
boundaries = [0, 117, 140, 183]       ## Converted from HbA1c boundaries: Normal is <5.7 %, prediabetes 5.7-<6.5, diabetes 6.5, poorly controlled diabetes >= 8.0%
labels_diabetes = ["", "Pre-diabetes", "Diabetes", "Poorly-controlled diabetes"]    ## Removed the Normal range label as it was implying 0 is normal range

# Draw dotted horizontal lines at the diabetes boundaries
for i, y in enumerate(boundaries[1:], start=1):                                             ## Iterate over the boundaries starting from 1 to avoid 0
    ax.axhline(y, color=crm.grayC(0.35), linestyle=(0,(2,2)), linewidth=1.2, alpha=1.0,     ## Draw dotted horizontal line at each boundary
               label="General blood glucose monitoring ranges" if i == 1 else "_nolegend_")   ## Add label of first line only to legend, to remove need of full label on y axis

# Add a right hand y-axis showing the diabetes thresholds at the band midpoints
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())  # Keeps both axes aligned

current_top = ax.get_ylim()[1]
band_tops = boundaries[1:] + [current_top]
mids = [(lo + hi) / 2.0 for lo, hi in zip(boundaries, band_tops)]

ax2.set_yticks(mids)                            ## Set right y-axis ticks to middle of points
ax2.set_yticklabels(labels_diabetes)            ## Labels for the right y-axis ticks
ax2.tick_params(axis="y", length=0, pad=14)     ## Ticks are 0 length, padding for labels to shift further right for manual brackets
#ax2.set_ylabel("Diabetes-related thresholds")   ## Right y-axis title, though there's the legend so not needed
ax2.grid(False)

## Add single artificial illustrative rounding example
## Use a fixed HbA1c to compute a realistic ±0.05% rounding effect in eAG.
example_hba1c_for_scale = 7.0
eag_low_ex = ngsp_to_eag_formula(example_hba1c_for_scale - 0.05)   ## eAG if example was rounded lower
eag_high_ex = ngsp_to_eag_formula(example_hba1c_for_scale + 0.05)   ## eAG if example was rounded higher
example_half_height = (eag_high_ex - eag_low_ex) / 2.0              ## Half vertical span giving symmetric +- error magnitude for the error bar

## Place the example near the legend - REMEMBER only adjust this based on saved image, not the interactive one as that renders live and isn't where the example actually is.
x_example = 4.61
y_example = 255 ## 240 for if the 90% PI bar is present

## Draw the example point estimate marker with a blue vertical ± error bar.
ax.errorbar(
    x_example, y_example,                                       ## Location of the point using the previous positioners
    yerr=example_half_height,                                   ## Symmetrical ± y dimension error bar
    fmt="none", # markersize=0.5, color="black",                ## Marker style, size and colour
    ecolor=crm.bamako(0.62), elinewidth=0.6,                    ## Error bar colour and thickness
    capsize=3, capthick=1.0,                                    ## Add caps to try and make it clearer it's an error bar
    alpha=1.0, zorder=6,                                        ## Draw on top of other bits (might need adjusting if other changes)
    label="_nolegend_"                                          ## Label for legend Example uncertainty due to rounding (±0.05 HbA1c)
)

## Add marker label same size and color as legend text
ax.text(
    4.7, 255,   ## 240 for if the 90% PI bar is present
    "Rounding uncertainty ± 0.05 HbA1c (to scale)",
    fontsize=8, color="black",
    ha="left", va="center"
)

## Axis titles, labels, title and legend
ax.set_xlabel("HbA1c (%)")
ax.set_ylabel("eAG (mg/dL)")
ax.set_title("Contributors to uncertainty when converting HbA1c to estimated average glucose")
leg = ax.legend(                ## Fix legend to lower right hand corner to permit dialling in the example rounding error point and label
    loc="upper left",
    frameon=False,
    handlelength=1.5,
    handletextpad=0.6,
    fontsize=8
)
#plt.suptitle("Data show conversion model uncertainty only; assay jitter & rounding are overlays", fontsize=10, y=0.98)
plt.annotate(f"Source: {filename}", xy=(0.01, 0.01), xycoords="figure fraction", fontsize=9, color="gray")  ## Put plot's data source filename in bottom left corner for traceability

output_filename = filename.replace(".csv", "_violin_plot_Figure_1.png")  ## Use same name as source file but swap .csv for plot name, for source traceability

plt.tight_layout()                                  ## Adjust to prevent labels being cut off
plt.savefig(output_filename, dpi=400)         ## Save at decent res for future zooming
print(f"Plot saved as: {output_filename}")          ## Console confirmation of save
plt.show()
