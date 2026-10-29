"""
Description:
Provides the calculations for the 90 % prediction intervals and the t-distribution samples.

Glossary:
eAG/eag: Estimated Average Glucose
HbA1c/hba1c: Glycated haemoglobin/hemoglobin
NGSP: National Glycohemoglobin Standardization Program

Author: Hugh Claridge
Last edited date: 29th October 2026
Version: v1.0
Version history: Public version of the statistics tools module
"""

import numpy as np
from scipy.stats import t

## Generate t-distributed samples
def generate_t_distribution_samples(hba1c_values, formula, beta_1=4.81, beta_2=2.03, num_samples=1000, n=507, rng=None):
    """
    Generates t-distributed samples for each HbA1c-based eAG value using Nathan et al.'s variance formula.
    :param hba1c_values: List of HbA1c values - which type depends on the formula which is next but is % for this work, but can be changed.
    :param formula: The conversion formula (i.e. ngsp_to_eag_formula for this work, but can be changed).
    :param beta_1: Variance scaling constant which uniformly adjusts the variance scaling (magnitude of the variance, up/down).
    :param beta_2: Variance growth exponent controlling the rate at which the variance grows with growing HbA1c.
    :param num_samples: Default number of samples to generate per eAG value.
    :param n: Sample size of 507 from Nathan et al. 2008's study size, it is not something we change with how many we are simulating as it is about the original uncertainty Nathan found and is part of the scale (SD) of t-distribution.
    :param rng: Use call's random number generator if specified, but if not, make an rng (rng always called in this project for reproducibility).
    :return: List of lists - one inner list per input HbA1c, each containing X t-distributed eAG samples centered on that row’s eAG point estimate (dictionaries bad as overwrite)
    :explainer: Was previously a dictionary but changed to list of lists because as a dictionary, eAG was a key, and dicts. only allow unique keys, so it would overwrite others
    :explainer: This explains why for multiple instances of e.g. 5.3 and 5.3, T1-T10 were identical.
    """
    rng = rng or np.random.default_rng()    ## If RNG supplied keep it, but if left as None, make one
    df = num_samples - 1    ## Degrees of freedom for t-distribution based on sample count
    row_t_samples = []      ## List collecting list of t-values per HbA1c

    for hba1c in hba1c_values:  ## Loop through each input HbA1c
        eag = formula(hba1c)    ## Point estimate eAG for the input HbA1c, centre of the t-distribution

        ## Compute scale dynamically using Nathan et al.'s formula: √(1+(1/n)) · √β₁ · HbA1c^(β₂/2)
        scale = np.sqrt(1 + (1 / n)) * np.sqrt(beta_1) * (hba1c ** (beta_2 / 2))

        ## Draw num_samples from the Student t-distribution centred on the point estimate (eag) with scale=horizontal stretch factor
        samples = eag + scale * rng.standard_t(df, size=num_samples)     ## Array of the defined number of t-distribution samples for the eAG
        row_t_samples.append(samples.tolist())                          ## Convert array to a list and append to the list of row t samples

    return row_t_samples      ## List of lists - for each HBA1c (outer list), an inner list of their t-dist draws

## Compute 90 % Prediction Interval
def compute_prediction_interval(hba1c_values, formula, beta_1=4.81, beta_2=2.03, n=507, alpha=0.1):
    """
    Generates the 90 % prediction interval for eAG given HbA1c, using Nathan et al. (2008).
    Paper link: https://doi.org/10.2337/dc08-0545
    Don't forget that alpha and beta are used a few times in the publication and not for the same values!
    :param hba1c_values: List of HbA1c values - which type depends on the formula which is next but it's % for this project, but can be changed for others in core_formulae
    :param formula: The conversion formula (i.e. ngsp_to_eag_formula for this project, but can be changed for others in core_formulae).
    :param beta_1: Variance scaling constant which uniformly adjusts the variance scaling (magnitude of the variance, up/down).
    :param beta_2: Variance growth exponent controlling the rate at which the variance grows with increasing HbA1c.
    :param n: Sample size (default here is 507 from Nathan et al. 2008).
    :param alpha: Significance level for the prediction interval (0.1 = 90 % CI, as tail probability per side is alpha/2 = 0.05).
    :return: Two parallel/paired arrays of the lower and upper bounds for the prediction interval per HbA1c provided.
    :explainer: t_critical is the critical t-value for a 90 % two-tailed interval.
    :explainer: The cut-off multiplied by the Standard Error (SE) to build a confidence or prediction interval is called the 'critical value'.
    :explainer: As this is a Student t-dist it's t-critical. This gives the half width.
    :explainer: t.ppf is the percent-point function (inverse of the cumulative distribution function CDF (?)) in the 't' object from scipy-stats
    """
    t_critical = t.ppf(1 - alpha / 2, df=n - 1)     ## Two-tailed t-score for 90 % prediction interval
    lowers, uppers = [], []                         ## Parallel arrays of lower and upper bounds for each Hba1c

    for hba1c in hba1c_values:          ## Loop through each HbA1c result
        eag = formula(hba1c)            ## Point estimate eAG for the HbA1c result

        ## Compute prediction SD using Nathan et al.'s formula: √(1+(1/n)) · √β₁ · HbA1c^(β₂/2)
        sd_pred = np.sqrt(1 + (1 / n)) * np.sqrt(beta_1) * (hba1c ** (beta_2 / 2))

        ## Multiply prediction SD by t-critical to get half-width and subtract/add from/to eAG point estimate to form 90 % prediction interval
        lowers.append(eag - t_critical * sd_pred)              ## Bottom edge of 90 % PI
        uppers.append(eag + t_critical * sd_pred)              ## Upper edge of 90 % PI

    return lowers, uppers     ## Parallel/paired arrays of lower/upper bounds
