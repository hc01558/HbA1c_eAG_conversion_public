"""
Description:
Contains useful functions for use throughout the conversion project (not all used, but retained for future use).

Author: Hugh Claridge
Start date: January 2025
Last edited date: 29th October 2026
Version: v1.0
Version history: Public version of the converter tools module
"""
import numpy as np
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP
from datetime import datetime

def generate_data(distribution_func, num_values=None, **kwargs):
    """
    Generates input data using the input distribution.

    :param distribution_func: Distribution function to generate data with different distributions e.g. uniform, gaussian etc.
    :param num_values: Number of values to generate. Added =None for the linear distribution's Range mode option.
    :param kwargs: Additional arguments for distribution generation (mean, std, etc.).
    :else: Added to be compatible with the fixed linear distribution's Range mode option which doesn't have a defined num_values.
    :return: Generated data.
    """
    if num_values is None:                                  ## If when called num_values left as None e.g. for linear dist's Range mode then...
        return distribution_func(**kwargs)                  ## Call the distribution function as is.
    else:                                                   ## Otherwise (i.e. because it wasn't left as None)
        return distribution_func(num_values, **kwargs)      ## Pass the specified inputs

def apply_percentage_jitter(values, percent=0.023, rng=None):
    """
    Applies percentage-based uncertainty (jitter) to input values i.e. as an input value increases so can the uncertainty
    e.g. 2.3% on 5 can go to 5.1 but on 10 can go to 10.2
    Got 0.023 as standard uncertainty, by halving the 0.046 expanded uncertainty from the +-4.6% reported by Unal and Atikeler 2018. This isn't used in the plots for the paper.
    Paper link: 10.14744/ijmb.2017.76486

    :param percent: Relative uncertainty (e.g., 0.023 for ±2.3%)
    :param rng: Use call's random number generator if specified, but if not, make an rng (rng always called in this project for reproducibility).
    :param values: List of input values (here HbA1c %)
    :return: Array of jittered values
    """
    rng = rng or np.random.default_rng()                         ## If RNG supplied keep it, but if left as None, make one
    noise = rng.normal(loc=0.0, scale=percent, size=len(values)) ## Ask the RNG to generate random numbers from a bell curve, centered on mean of 0.0, with spread (SD) a defined %, of the length one number for each element in values
    return values * (1 + noise)                                  ## Now multiply the defined value by its noise, thus noise scales with magnitude

def apply_sd_jitter(values, sd=0.15, rng=None):
    """
    Applies SD-based uncertainty (jitter) to input values, each being impacted by a random amount drawn from 0 sd +-0.15 i.e. doesn't increase as values increase.
    This isn't used in the plots for the paper.

    :param values: List of input values (here HbA1c %)
    :param sd: Standard deviation of uncertainty to add
    :param rng: Use call's random number generator if specified, but if not, make an rng (rng always called in this project for reproducibility).
    :return: Array of jittered values
    """
    rng = rng or np.random.default_rng()                        ## If RNG supplied keep it, but if left as None, make one
    noise = rng.normal(loc=0.0, scale=sd, size=len(values))     ## Ask the RNG to generate random numbers from a bell curve, centered on mean of 0.0, with spread (SD) a defined number, of the length one number for each element in values
    return values + noise                                       ## Now add the noise, so noise does not scale with magnitude

def convert_values(values, formula):
    """
    Applies the conversion formula to the given list of values.
    :param values: List of input values.
    :param formula: Formula function which applies the defined conversion formula to each input value.
    :return: List of converted values.
    """
    return [formula(val) for val in values]               ## Calls the defined formula function on each value, for each value in the list of values, and returns the outputs into a new list

def round_half_up(x, ndigits=0):
    """
    Lab behaviour of rounding a number to the nearest unit, with .5 rounding up (away from zero).
    This is because the standard inbuilt round() behaviour is banker's round to even (e.g. .45 rounds down to .4 whilst .55 rounds up to .6)

    :param x: The input value to round (float/int/Decimal).
    :param ndigits: Number of decimal places to round to: default is 0 = integer (e.g. eAG mg/dL), in calls often 1 = 1 dp (e.g. HbA1c %)
    :return: Rounded value as a float.
    """
    if pd.isna(x):                  ## NaN safe keeping missings as Nan
        return np.nan
    quantum = Decimal('1').scaleb(-ndigits)                     ## Builds the step/quantum to round to, here starting at 1.0, and then move dp -ndigits, so 0 means int, 1 means 0.1 etc.
    d = x if isinstance(x, Decimal) else Decimal(str(x))        ## If value is decimal then use as is, else convert to string (otherwise e.g. 6.45 as a float is often actually stored as 6.44999... thus could round wrong way if used directly)
    return float(d.quantize(quantum, rounding=ROUND_HALF_UP))   ## Rounds value d to the nearest defined step/quantum and convert decimal back to float for compatibility (default for Decimal is again banker's rounding to even)

def construct_filename(input_values, num_samples=None):
    """
    To remove the need to manually edit file names, this makes one based on number of input values, number of t-samples, and input range.
    :param input_values: The input values (almost always HbA1c % but could change).
    :param num_samples: The number of t-distribution samples per input
    :return: A string filename with timestamp and .csv extension.
    """
    num_values = len(input_values)                               ## Count how many results
    min_val = round_half_up(min(input_values), 1)         ## Find the lowest value, rounded to 1 decimal place using our rounding
    max_val = round_half_up(max(input_values), 1)         ## Find the highest value, rounded to 1 decimal place using our rounding
    sample_str = f"s{num_samples}" if num_samples is not None else "sNA"        ## If num_samples provided then write it, if blank then write sNA
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")                        ## Date+time when file made
    return f"{timestamp}_n{num_values}_{sample_str}_{min_val:.1f}-{max_val:.1f}.csv"   ## Join bits together with _ between

def save_to_csv(original_values, input_values, converted_values, lower_bounds, upper_bounds, t_samples,
                filename="converted_values.csv"):
    """
    Saves the original HbA1c, lab-style rounded input HbA1c, converted eAG point estimates, and simulation results to a .csv file.
    Checks that the columns are the same length, and then saves the defined columns. T-sample cols auto-adjust given how often this needed changing.

    :param original_values: Original HbA1c values before rounding or jitter
    :param input_values: Rounded and/or jittered HbA1c values actually used in conversion, given that's what someone will actually convert
    :param converted_values: Resulting eAG point estimate values, unrounded
    :param lower_bounds: Prediction interval lower bounds, unrounded
    :param upper_bounds: Prediction interval upper bounds, unrounded
    :param t_samples: List of lists containing t-distribution samples for each HbA1c row
    :param filename: Output CSV file name
    """
    ## Compares length of input_values to other cols to ensure they match its length and flag if something doesn't
    n = len(input_values)
    assert len(original_values) == n, "original_values length mismatch"
    assert len(converted_values) == n, "converted_values length mismatch"
    assert len(lower_bounds) == n, "lower_bounds length mismatch"
    assert len(upper_bounds) == n, "upper_bounds length mismatch"

    for i, col in enumerate(t_samples, 1):          ## Loop to check over the variable number of t sample cols
        assert len(col) == n, f"T-Sample-{i} length mismatch"

    data = {                                       ## Builds a dictionary with these keys the .csv col headers
        "Original HbA1c": original_values,
        "Input Values": input_values,
        "Converted Values": converted_values,
        "Lower Bound": lower_bounds,
        "Upper Bound": upper_bounds,
    }
    for i, col in enumerate(t_samples, 1):      ## Using a loop to create the name with t sample counter starting at 1...
        data[f"T-Sample-{i}"] = col             ## Add the column i.e. T-Sample-1, T-Sample-2 etc. making nicer column headers and auto adjustable given changing of number of t-draws

    pd.DataFrame(data).to_csv(filename, index=False)   ## Default is to have row index as first column, 0, 1, 2, etc. this disables that as it's not useful
    print(f"Data saved to {filename}")                 ## Console message about the save