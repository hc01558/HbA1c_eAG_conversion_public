"""
Description:
This is the main script of the project which generates and saves a dataset of:
1. A defined number (here 10,000) of randomly drawn or fixed linear distribution of NGSP HbA1c % values within defined range (could be a different input unit too)
2. Assay uncertainty can be included here, but hasn't been used here in the end.
3. Round HbA1c to 1dp lab-style for use by the chosen conversion model.
4. Convert HbA1c to eAG (or from/to other units) using the defined conversion model result for each of the HbA1c % - these create the point estimates.
5. Use Nathan et al.'s 2008 variance model (http://dx.doi.org/10.2337/dc08-0545) to quantify the prediction uncertainty around each eAG point estimate, and draw t-distributed Monte-Carlo samples of eAG to simulate that error.
6. Calculation of 90 % prediction intervals for those eAG results using the variance model from Nathan et al. 2008.
7. Saves all results to a .csv file

The only rounding happens to the HbA1c % which are used in the conversion, to be like the results reported from a lab, and someone then wants to convert them to eAG.
Others remain unrounded to enable separating sources of uncertainty (rounding, assay jitter, conversion uncertainty).

Glossary:
eAG: Estimated Average Glucose
HbA1c: Glycated haemoglobin/hemoglobin
NGSP: National Glycohemoglobin Standardization Program

Author: Hugh Claridge
Start date: January 2025
Last edited date: 29th October 2026
Version: v1.0
Version history: Public version of the main module
"""
## Imports from helper modules
from core_converter_tools import generate_data, convert_values, round_half_up, construct_filename, save_to_csv, apply_percentage_jitter, apply_sd_jitter
from core_distributions import uniform_distribution, gaussian_distribution, fixed_linear_distribution
from core_formulae import ngsp_to_eag_formula
from core_statistics import compute_prediction_interval, generate_t_distribution_samples

import numpy as np

def main():
    ## Toggle and set the global seed and create random number generator (rng) - this is for reproducible randomness
    use_seed = True
    global_seed = 42
    rng = np.random.default_rng(global_seed) if use_seed else np.random.default_rng()

    ## Define the number of inputs (HbA1c % here) to generate
    num_values = 10000

    ## Use this for uniform distribution (adjust low/high for desired range but remember if rounding, starting at e.g. 5.1 means you only get 5.1 to 5.15, then 5.15-5.24)
    hba1c_values = generate_data(uniform_distribution, num_values=num_values, low=4.45, high=9.45, rng=rng)
    ## Use this for gaussian distribution for ~99.7 % (+-3 SD) coverage of 5.1-5.7
    #hba1c_values = generate_data(gaussian_distribution, num_values=num_values, mean=5.4, std=0.1, rng=rng)
    ## Use this for gaussian distribution for ~99.7 % (+-3 SD) coverage of 4.5-8.0ish
    # hba1c_values = generate_data(gaussian_distribution, num_values=num_values, mean=6.25, std=0.583, rng=rng)
    ## Use this for gaussian distribution for ~99.7 % (+-3 SD) coverage of 4.0-13.9ish
    # hba1c_values = generate_data(gaussian_distribution, num_values=num_values, mean=8.95, std=1.65, rng=rng)
    ## Use this for fixed linear distribution
    # hba1c_values = generate_data(fixed_linear_distribution, num_values=None, low=4.0, high=14.0, step=0.1)

    ## Check HbA1c value generation
    assert len(hba1c_values) == num_values, "Input count doesn't equal num_values"
    assert np.isfinite(hba1c_values).all(), "NaNs or infs in inputs"                               ## IsFinite checks both
    assert hba1c_values.min() >= 4.45 and hba1c_values.max() < 9.45, "Raw inputs generated outside of 4.45 to 9.45"

    ## Store the original HbA1c values pre-rounding (and sometimes pre-jitter) for future comparison
    original_hba1c_values = hba1c_values.copy()

    ## Can apply assay uncertainty here to the generated HbA1c values or not, from Unal and Atikeler 2018
    ## Percent-based jitter scales with magnitude, SD is 2.3 % of the result (CV = 2.3 %)
    # hba1c_values = apply_percentage_jitter(hba1c_values, percent=0.023, rng=rng)

    ## Jitter is causing HbA1c % 'inputs' to fall outside of the defined range, so can clip here
    # hba1c_values = np.clip(hba1c_values, 4.0, 13.9)

    ## Now round to HbA1c %'s 1 decimal place, i.e. how test results are reported
    hba1c_values = [round_half_up(v, 1) for v in hba1c_values]  ## Lab-style round half up instead of inbuilt banker's rounding
    assert min(hba1c_values) >= 4.5 and max(hba1c_values) <= 9.4, "Rounded inputs outside of 4.5 to 9.4"

    ## Convert HbA1c values to eAG point estimates
    converted_values = convert_values(hba1c_values, ngsp_to_eag_formula)

    ## Check converted values generation
    assert len(converted_values) == num_values, "Count of converted values doesn't match inputs"
    assert np.isfinite(converted_values).all(), "NaNs or infs in converted values"

    ## Generate t-distribution samples for each eAG value
    ## Define the number of t-distribution samples
    num_samples = 1000
    row_t_samples = generate_t_distribution_samples(
        hba1c_values, formula=ngsp_to_eag_formula, num_samples=num_samples, rng=rng
    )

    ## Transpose t-samples from per row lists i.e. each HbA1c's draws, to per t-sample nth draw cols
    t_samples_columns = list(map(list, zip(*row_t_samples)))

    ## Check t-sample generation
    # Jess noted 'assert commands that are in loops show errors', but I suspect it was a linter warning as I cannot reproduce it.
    assert len(t_samples_columns) == num_samples, "Col count of t-samples doesn't equal num_samples"
    for i, col in enumerate(t_samples_columns, 1):
        arr = np.asarray(col)
        assert arr.shape == (num_values,), f"T-Sample-{i} col length doesn't equal num_values"
        assert np.isfinite(arr).all(), f"NaNs or infs in T-Sample-{i} col"

    ## Calculate 90 % prediction intervals (paired/parallel arrays of lower and upper bounds for each HbA1C)
    lower_bounds, upper_bounds = compute_prediction_interval(hba1c_values, formula=ngsp_to_eag_formula)

    ## Check prediction intervals
    assert np.isfinite(lower_bounds).all(), "NaNs or infs in lower PI bounds"
    assert np.isfinite(upper_bounds).all(), "NaNs or infs in upper PI bounds"
    assert (lower_bounds < upper_bounds).all(), "Lower > upper in prediction interval"

    ## Check to see if the simulated t-sample 90 % coverage intervals reproduce the analytic 90 % prediction intervals
    t_matrix = np.column_stack(t_samples_columns)  # shape num_values, num_samples
    sim_p5 = np.quantile(t_matrix, 0.05, axis=1)
    sim_p95 = np.quantile(t_matrix, 0.95, axis=1)

    max_abs_low = float(np.max(np.abs(sim_p5 - lower_bounds)))
    max_abs_high = float(np.max(np.abs(sim_p95 - upper_bounds)))

    print(
        f"Checking t-samples vs analytic 90 % PI: "
        f"max |P5 - Lower| = {max_abs_low:.3f} mg/dL; "
        f"max |P95 - Upper| = {max_abs_high:.3f} mg/dL"
    )

    abs_low = np.abs(sim_p5 - lower_bounds)
    abs_high = np.abs(sim_p95 - upper_bounds)

    print(
        "Abs diff summary (P5 vs Lower): "
        f"median={np.median(abs_low):.3f}, "
        f"p95={np.percentile(abs_low, 95):.3f}, "
        f"max={np.max(abs_low):.3f} mg/dL"
    )
    print(
        "Abs diff summary (P95 vs Upper): "
        f"median={np.median(abs_high):.3f}, "
        f"p95={np.percentile(abs_high, 95):.3f}, "
        f"max={np.max(abs_high):.3f} mg/dL"
    )

    ##
    #Checking t-samples vs analytic 90% PI: max |P5 - Lower| = 4.624 mg/dL; max |P95 - Upper| = 5.200 mg/dL
    #Abs diff summary (P5 vs Lower): median=0.707, p95=2.147, max=4.624 mg/dL
    #Abs diff summary (P95 vs Upper): median=0.688, p95=2.138, max=5.200 mg/dL

    ## Save everything using the construct_filename and save_to_csv functions - decimal point rounding control of the converted values is in converter_tools.py
    save_to_csv(
        original_values=original_hba1c_values,
        input_values=hba1c_values,
        converted_values=converted_values,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        t_samples=t_samples_columns,
        filename=construct_filename(hba1c_values, num_samples=num_samples)
    )

if __name__ == "__main__":          ## Only run this module when run directly
    main()
