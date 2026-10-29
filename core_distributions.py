"""
Description:
Contains the various distribution functions that can be used when generating data.
We use uniform for the project, but explored Gaussian and linear too; functionality retained in case useful in future, but other code will need changing to work with these (sometimes available commented-out where previously used).

Author: Hugh Claridge
Start date: January 2025
Last edited date: 11th November 2025
Version: v1.0
Version history: Added lines for author, last edit date, version and version history as per Jess's review suggestion.
"""
import numpy as np

## Uniform distribution
def uniform_distribution(num_values, low=0, high=15, rng=None):
    """
    Draws values from a uniform distribution - can look non-uniform with lower numbers, so sense check with 10,000+ in this project
    :param num_values: Number of values to generate, has to be provided
    :param low: Lowest value in the range, default provided
    :param high: Highest value in the range, default provided
    :param rng: Use call's random number generator if specified, but if not, make an rng (rng always called in this project for reproducibility)
    :return: An array of length num_values between specified or default lower and upper bounds with equal probability
    """
    rng = rng or np.random.default_rng()
    return rng.uniform(low, high, num_values)

## Gaussian/normal distribution
def gaussian_distribution(num_values, mean=0, std=1, rng=None):
    """
    Draws values from a gaussian/normal distribution - explored this but got messy so reverted to uniform as clearer to communicate
    :param num_values: Number of values to generate, must be provided
    :param mean: The distribution mean, output centred around this.
    :param std: Distribution standard deviation, controls spread.
    :param rng: Use call's random number generator if specified, but if not, make an rng (rng always called in this project for reproducibility)
    :return: An array of length num_values with bell-curve distribution centered on specified or default mean and SD
    """
    rng = rng or np.random.default_rng()
    return rng.normal(mean, std, num_values)


def fixed_linear_distribution(num_values=None, low=4.0, high=14.0, step=0.1):
    """
    Draws values from a fixed linear distribution - useful for early exploratory analysis with small numbers to avoid distribution-created issues
    Range mode and Count modes for different generation requirements.
    :param num_values: Number of values to generate, must be provided - just so happens to be 100 for this range of 4.0 - 14 at 0.1 increments
    :param low: Lowest value in the range, default provided
    :param high: Highest value in the range, default provided
    :param step: Defines the increment
    :return: An array of values to 1dp (HbA1c % lab report style inputs)
    """
    if num_values is None:
        ## Range mode: num_values not provided, fixed step, length determined by low and high
        values = np.arange(low, high, step)
    else:
        ## Count mode: num_values defines length, low defines starting point, high ignored, increment determined by step
        values = low + step * np.arange(num_values)
    return values


