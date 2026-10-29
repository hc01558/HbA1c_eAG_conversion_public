"""
Description:
Contains various blood glucose and HbA1c point estimate/deterministic conversion formulae which can be used in the project.
Only used NGSP/DCCT HbA1c (%) to eAG for this project but others may be useful in future thus retained.

Glossary:
eAG: Estimated Average Glucose
HbA1c: Glycated haemoglobin/hemoglobin
IFCC: International Federation of Clinical Chemistry
NGSP: National Glycohemoglobin Standardization Program

Author: Hugh Claridge
Start date: January 2025
Last edited date: 29th October 2026
Version: v1.0
Version history: Public version of the formulae module
"""

## NGSP HbA1c (%) to eAG (mg/dL) conversion formula
def ngsp_to_eag_formula(x):
    """
    Converts NGSP HbA1c (%) to eAG (mg/dL) using the formula from Nathan et al. 2008
    Paper link: http://dx.doi.org/10.2337/dc08-0545
    :param x: The input NGSP HbA1c (%) value
    :return: The output eAG (mg/dL) value
    """
    a = 28.7
    b = 46.7
    return (a * x) - b

## eAG (mg/dL) to NGSP HbA1c (%) conversion formula
def eag_to_ngsp_formula(x):
    """
    Converts eAG (mg/dL) to NGSP HbA1c (%) using the reverse of the NGSP to eAG formula from Nathan et al. 2008
    Paper link: http://dx.doi.org/10.2337/dc08-0545
    :param x: The input eAG (mg/dL) value
    :return: The output NGSP HbA1c (%) value
    """
    a = 28.7
    b = 46.7
    return (x + b) / a

## IFCC HbA1c (mmol/mol) to NGSP HbA1c (%) conversion formula
def ifcc_to_ngsp_formula(x):
    """
    Converts IFCC HbA1c (mmol/mol) to NGSP HbA1c (%) using the formula from IFCC.
    Paper link: https://ngsp.org/docs/IFCCstd.pdf
    :param x: The input IFCC HbA1c (mmol/mol) value
    :return: The output NGSP (%) value
    """
    a = 0.09148
    b = 2.152
    return (a * x) + b

## NGSP HbA1c (%) to IFCC HbA1c (mmol/mol) conversion formula
def ngsp_to_ifcc_formula(x):
    """
    Converts NGSP HbA1c (%) to IFCC HbA1c (mmol/mol) using the formula from IFCC.
    Paper link: https://ngsp.org/docs/IFCCstd.pdf
    :param x: The input NGSP HbA1c (%) value
    :return: The output IFCC HbA1c (mmol/mol) value
    """
    a = 10.93
    b = 23.50
    return (a * x) - b

## NGSP HbA1c (%) to IFCC HbA1c (%) conversion formula
def ngsp_to_ifcc_percent_formula(x):
    """
    Converts NGSP HbA1c (%) to IFCC HbA1c (%) using the formula from Sacks et al. 2005.
    This is extremely uncommon now given the confusion it causes, with IFCC presentation standardised as mmol/mol long ago, for good reason.
    Sacks et al. published lower significant figures for this formula vs. the IFCC HbA1c (mmol/mol) to NGSP HbA1c (%), so have used the formula's more granular numbers than Sacks' paper.
    Paper link: https://doi.org/10.1373/clinchem.2004.047431
    :param x: The input NGSP HbA1c (%) value
    :return: The output IFCC HbA1c (%) value
    """
    a = 0.9148
    b = 0.02152
    return (a * x) + b
