"""
Description:
Calculates and saves a .csv of the per 0.1 % HbA1c bin misclassification summaries for point estimates and t-draws in same/different classification bands.

Glossary:
eAG: Estimated Average Glucose
HbA1c: Glycated haemoglobin/hemoglobin
NGSP: National Glycohemoglobin Standardization Program

Author: Hugh Claridge
Start date: January 2025
Last edited date: 29th October 2026
Version: v1.0
Version history: Public version of the misclassification likelihood analysis
"""

import pandas as pd
import numpy as np
from core_converter_tools import round_half_up

## Load in the spreadsheet for checking and track the filename for plot tracing
## By default this will run on the data file used in the paper
filename = "20260104_233659_n10000_s1000_4.5-9.4.csv" ## Corrected scale formula and 90% PI formula and rng for t-dist samples

## Load chosen CSV file into dataframe
df = pd.read_csv(filename)

## Create list of T-sample column names
t_cols = [c for c in df.columns if c.startswith("T-Sample-")]

## Using the list of t_sample cols, calculate the max so the final bin can be closed (added to resolve issue when dealing with ranges converting to above 161 mg/dL)
max_val = df[t_cols].max().max()    ## Take the t-sample cols, first max is column-wise highest eAG, second max is highest of that series, so highest of the highest

## Binning with the right and left boundaries for diabetes mg/dL boundaries
last_edge = max(max_val + 1, 162)   ## This is to permit using varying ranges which can convert to lower or higher eAGs beyond the diabetes boundaries
bins = [0, 117, 140, 183, last_edge]  ## Converted from HbA1c boundaries: Normal is <5.7 %, prediabetes 5.7-<6.5, diabetes 6.5, poorly controlled diabetes >= 8.0%
labels = ["Normal", "Pre-diabetes", "Diabetes", "Poorly-controlled diabetes"] ## Diabetes categories by eAG range - potential implication that 0 is seen as 'normal'

## Function to assign values (converted/point estimates and t-samples) to diabetes categories
## pd.cut takes the numeric value column and assigns it to the defined categorical bin, thus creates a categorical series where each row contains the diabetes category of the value (whether point estimate or t-sample)
## pd.cut is a direct comparison so 99.99999 remains in the lower bin and is not taken into 100. Intentional non-rounding here as we are only looking at conversion model issues rather than introducing rounding.
def assign_to_band(col):
    return pd.cut(col, bins=bins, labels=labels, right=False, ordered=True) ## ordered=True preserves the left to right order of labels so that it goes normal, pre-diabetes, diabetes etc. Needed to permit measuring how many categories away from point estimate

## Assign each row's point-estimate (Converted Value) into a diabetes category
## This is the 'reference' band t-samples will be compared against
reference_band = assign_to_band(df["Converted Values"])

## Assign each row's t-sample values into diabetes categories
banded_t_samples = df[t_cols].apply(assign_to_band)

## Calculate the % of the row’s t-samples that land in the same category as the point estimate (distance = 0), 1 category away (distance =1) etc.
## Need to convert the ordered categorical labels to integer codes for subtraction
## Code mapping given the labels of Normal = 0, Pre-diabetes = 1 etc.
category_codes = reference_band.cat.codes
sample_category_codes = banded_t_samples.apply(lambda c: c.cat.codes)

## Calculate category distance for every t-sample in every row (for separate category +1 above, -1 below etc.)
category_distance = sample_category_codes.sub(category_codes, axis=0) ## Creates data frame with one row per input and one column per t-sample draw
## Calculate absolute distance for every t-sample in every row (for combined series ±1 categories, ±2 etc. and total misclassification
abs_distance = category_distance.abs()

## Convert distances into row level percentages (0–100) for each distance bucket
## Separate directional distances for separate line plots
percent_same_category = (category_distance == 0).mean(axis=1) * 100.0
percent_one_category_above = (category_distance == 1).mean(axis=1) * 100.0
percent_one_category_below = (category_distance == -1).mean(axis=1) * 100.0
percent_one_or_more_category_above = (category_distance >= 1).mean(axis=1) * 100.0
percent_one_or_more_category_below = (category_distance <= -1).mean(axis=1) * 100.0
## Separate directional distances 2 or more
percent_two_category_above = (category_distance == 2).mean(axis=1) * 100.0
percent_two_category_below = (category_distance == -2).mean(axis=1) * 100.0
percent_two_or_more_category_above = (category_distance >= 2).mean(axis=1) * 100.0
percent_two_or_more_category_below = (category_distance <= -2).mean(axis=1) * 100.0
percent_three_or_more_category_above = (category_distance >= 3).mean(axis=1) * 100.0
percent_three_or_more_category_below = (category_distance <= -3).mean(axis=1) * 100.0

## Combined (absolute) distances
percent_abs_1 = (abs_distance == 1).mean(axis=1) * 100.0
percent_abs_1plus = (abs_distance >= 1).mean(axis=1) * 100.0       ## This is the same as overall misclassificaiton
percent_abs_2 = (abs_distance == 2).mean(axis=1) * 100.0
percent_abs_2plus = (abs_distance >= 2).mean(axis=1) * 100.0
percent_abs_3plus = (abs_distance >= 3).mean(axis=1) * 100.0

## Overall missclassification
percent_misclass = (abs_distance >= 1).mean(axis=1) * 100.0

## Store both raw and rounded versions for directional category distances
## Raw
df["Same % (raw)"] = percent_same_category
df["+1 % (raw)"] = percent_one_category_above
df["-1 % (raw)"] = percent_one_category_below
df["≥+1 % (raw)"] = percent_one_or_more_category_above
df["≤-1 % (raw)"] = percent_one_or_more_category_below
df["+2 % (raw)"] = percent_two_category_above
df["-2 % (raw)"] = percent_two_category_below
df["≥+2 % (raw)"] = percent_two_or_more_category_above
df["≤-2 % (raw)"] = percent_two_or_more_category_below
df["≥+3 % (raw)"] = percent_three_or_more_category_above
df["≤-3 % (raw)"] = percent_three_or_more_category_below
## Rounded
df["Same %"] = df["Same % (raw)"].apply(lambda x: round_half_up(x, 1))
df["+1 %"] = df["+1 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["-1 %"] = df["-1 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≥+1 %"] = df["≥+1 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≤-1 %"] = df["≤-1 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["+2 %"] = df["+2 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["-2 %"] = df["-2 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≥+2 %"] = df["≥+2 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≤-2 %"] = df["≤-2 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≥+3 %"] = df["≥+3 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["≤-3 %"] = df["≤-3 % (raw)"].apply(lambda x: round_half_up(x, 1))

## Store both raw and rounded versions for combined category distances
## Raw
df["±1 % (raw)"] = percent_abs_1
df["±1+ % (raw)"] = percent_abs_1plus          ## This is the same as overall misclassificaiton
df["±2 % (raw)"] = percent_abs_2
df["±2+ % (raw)"] = percent_abs_2plus
df["±3+ % (raw)"] = percent_abs_3plus
df["Misclass % (raw)"] = percent_misclass
## Rounded
df["±1 %"] = df["±1 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["±1+ %"] = df["±1+ % (raw)"].apply(lambda x: round_half_up(x, 1))
df["±2 %"] = df["±2 % (raw)"].apply(lambda x: round_half_up(x, 1))
df["±2+ %"] = df["±2+ % (raw)"].apply(lambda x: round_half_up(x, 1))
df["±3+ %"] = df["±3+ % (raw)"].apply(lambda x: round_half_up(x, 1))
df["Misclass %"] = df["Misclass % (raw)"].apply(lambda x: round_half_up(x, 1))

# Print rounded cols for checking
cols_to_show = [
    "Input Values", "Converted Values", "Same %", "+1 %", "-1 %", "≥+1 %", "≤-1 %", "+2 %", "-2 %", "≥+2 %", "≤-2 %","≥+3 %", "≤-3 %", "±1 %", "±1+ %", "±2 %","±2+ %", "±3+ %", "Misclass %"
]

print(df[cols_to_show].head(30).to_string(index=False))

## Aggregate by HbA1c input value (0.1 % bins)
bin_edges = np.arange(4.5, 9.6 + 1e-9, 0.1)            ## HbA1c bins: 4.5 - 9.50 at 0.1 increments (1e- needed due to binary floats messing with arange 9.0 edge)
df["HbA1c Bin"] = pd.cut(df["Input Values"], bins=bin_edges, right=False, include_lowest=True)        ## Assign each HbA1c input to a 0.1 % bin

## Per bin summary stats for the separate (i.e. not +-1 lumped together) bands
def _agg_mean(s): return s.mean()
def _agg_sd(s): return s.std(ddof=-1)  ## Sample SD rather than Population SD (ddof=0)

## Aggregation per 0.1 % HbA1c bin
## Mean columns are mean % of t-samples in each distance bucket (across rows in the HbA1c bin)
## SD columns are SD across rows of those percentages for the ribbons on the plot
## Aggregation is of raw columns only, the rounded cols only been for printing
bin_stats = df.groupby("HbA1c Bin", observed=True).agg(         ## Retains only categorical bins with data; false retains even empties
    Total=("Input Values", "size"),

    ## Separate directional distances
    MeanSame=("Same % (raw)", _agg_mean),
    SDSame=("Same % (raw)", _agg_sd),
    Mean1Above=("+1 % (raw)", _agg_mean),
    SD1Above=("+1 % (raw)", _agg_sd),
    Mean1Below=("-1 % (raw)", _agg_mean),
    SD1Below=("-1 % (raw)", _agg_sd),
    MeanMoreThan1Above=("≥+1 % (raw)", _agg_mean),
    SDMoreThan1Above=("≥+1 % (raw)", _agg_sd),
    MeanMoreThan1Below=("≤-1 % (raw)", _agg_mean),
    SDMoreThan1Below=("≤-1 % (raw)", _agg_sd),

    ## Two or more
    Mean2Above=("+2 % (raw)", _agg_mean),
    SD2Above=("+2 % (raw)", _agg_sd),
    Mean2Below=("-2 % (raw)", _agg_mean),
    SD2Below=("-2 % (raw)", _agg_sd),
    MeanMoreThan2Above=("≥+2 % (raw)", _agg_mean),
    SDMoreThan2Above=("≥+2 % (raw)", _agg_sd),
    MeanMoreThan2Below=("≤-2 % (raw)", _agg_mean),
    SDMoreThan2Below=("≤-2 % (raw)", _agg_sd),
    MeanMoreThan3Above=("≥+3 % (raw)", _agg_mean),
    SDMoreThan3Above=("≥+3 % (raw)", _agg_sd),
    MeanMoreThan3Below=("≤-3 % (raw)", _agg_mean),
    SDMoreThan3Below=("≤-3 % (raw)", _agg_sd),

    ## Combined (absolute) distances
    MeanAbs1=("±1 % (raw)", _agg_mean),
    SDAbs1=("±1 % (raw)", _agg_sd),
    MeanAbs1plus=("±1+ % (raw)", _agg_mean),        ## This is the same as overall mean misclassification
    SDAbs1plus=("±1+ % (raw)", _agg_sd),            ## This is the same as overall SD misclassification
    MeanAbs2=("±2 % (raw)", _agg_mean),
    SDAbs2=("±2 % (raw)", _agg_sd),
    MeanAbs2plus=("±2+ % (raw)", _agg_mean),
    SDAbs2plus=("±2+ % (raw)", _agg_sd),
    MeanAbs3plus=("±3+ % (raw)", _agg_mean),
    SDAbs3plus=("±3+ % (raw)", _agg_sd),

    ## Overall misclassification (1 or more category away)
    MeanMisclass=("Misclass % (raw)", _agg_mean),
    SDMisclass=("Misclass % (raw)", _agg_sd),
)

## Reset HbA1c bin from groupby's index to column for saving
bin_stats = bin_stats.reset_index()

## Check to see that's what's produced fits with expectations on the range
assert df["Input Values"].min() >= 4.5, "Minimum input below 4.5"
assert df["Input Values"].max() <= 9.4, "Maximum input above 9.4"
assert bin_stats["HbA1c Bin"].min().left >= 4.5, "Left bin edge below 4.5"
assert bin_stats["HbA1c Bin"].max().right >= 9.5, "Right bin edge below 9.5"

## Save the file
output_filename = filename.replace(".csv", "_misclass_analysis.csv")
bin_stats.to_csv(output_filename, index=False)
print(f"Analysis misclassification likelihood file saved as: {output_filename}")

## Show the rounded contents for quick check
cols = [c for c in bin_stats.columns if c.startswith(("Mean", "SD"))]
print_bin = bin_stats.copy()
print_bin[cols] = print_bin[cols].map(lambda x: round_half_up(x, 1))
print(print_bin.to_string(index=False))