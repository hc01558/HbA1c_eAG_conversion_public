## Quantifying prediction uncertainty and categorical misclassification risk in regression-based biomedical measurement conversions: an HbA1c to eAG case study
- This Python project produces the data files and plots for the published manuscript
- It currently generates glycated haemoglobin (NGSP HbA1c in %) inputs which are a measure of blood sugar, converts them to estimated average glucose (eAG in mg/dL) using Nathan et al.'s (2008) variance model, and quantifies prediction uncertainty with t-distributed samples and 90 % prediction intervals.
- Components can also be used to generate and convert to/from other blood glucose values, but do note that doing so will require more edits than simply changing to one of the other conversion formulae provided.
- It is in a modular structure to improve reusability.

Author: Hugh Claridge<br>
Start date: January 2025<br>
Last edited date: 29th October 2026<br>
Version: v1.0<br>
Version history: Public version of the code for the HbA1c to eAG conversion manuscript

Reference:<br>
Nathan et al. (2008): http://dx.doi.org/10.2337/dc08-0545

## Download 
You can clone this project, or use the below download link to obtain the source code zip.<br>
If you want to reproduce the plots as in the manuscript, you will also need to download the main simulation dataset and misclassification analysis dataset.<br>
Remember to save these data files into the extracted project folder, leaving the file names intact.<br>
Download link: https://github.com/hc01558/HbA1c_eAG_conversion_public/releases

## Glossary
- eAG: Estimated Average Glucose
- HbA1c: Glycated haemoglobin/hemoglobin
- IFCC: International Federation of Clinical Chemistry
- NGSP: National Glycohemoglobin Standardization Program

## Contents
- `core_main.py` — Main pipeline (generates data, performs conversion, produces t-samples and 90 % prediction intervals, saves a .csv).
- `core_converter_tools.py` — Helpers for the main pipeline (e.g. rounding, assay jitter, save, filename builder).
- `core_formulae.py` — Conversion formulae (NGSP HbA1c (%), IFCC HbA1c (mmol/mol and %), eAG (mg/dL) conversions).
- `core_statistics.py` — 90 % prediction intervals and generation of t-samples.
- `core_distributions.py` — Uniform/Gaussian/fixed-linear distributions.
- `analysis_misclassification_likelihood.py` — Per-bin misclassification analysis (in the same or a different diabetes category as the point estimate)
- `plot_violins_Figure_1.py` - Figure 1
- `plot_stepped_line_Figure_2_upper_panel.py` - Figure 2 upper panel
- `python plot_stepped_line_Figure_2_lower_panel.py` - Figure 2 lower panel

## Quick start
```python
# Run this to create a virtual environment to contain this project (not mandatory)
# Windows:
python -m venv .venv && .\.venv\Scripts\Activate.ps1

# Linux/macOS
python -m venv .venv && source .venv/bin/activate

# Run this to install required dependencies
pip install pandas numpy scipy matplotlib cmcrameri

# Create Figure 1 violin plot, using the study's own dataset
# Ensure you have downloaded the main dataset '20260104_233659_n10000_s1000_4.5-9.4.csv' and saved it into your project folder
python plot_violins_Figure_1.py

# Create Figure 2's two stepped line panels, using the study's own dataset
# Ensure you have downloaded the misclassification analysis data file '20260104_233659_n10000_s1000_4.5-9.4_misclass_analysis.csv' and saved it into your project folder
python plot_stepped_line_Figure_2_upper_panel.py
python plot_stepped_line_Figure_2_lower_panel.py

# If you would like to recreate the data files yourself: 
# Run the main data generation
python core_main.py

# Run misclassification analysis on the output .csv saved by core_main.py
python analysis_misclassification_likelihood.py

# You can then reuse the above plotting code for Figures 1 and 2
# However, you will need to replace the file names on which the plotting code runs with the new files' names

