# Replication code for Matching Techniques in Energy Policy Analysis: A Teaching Case for Smart Cities

This is the code repository for the Albany Teaching Case. The R-code file consists of the fixed-effects model and matching algorithms to compare energy use in Albany, GA, among HUD-funded programs' participants and non-participants. The data is deposited at [https://doi.org/10.7910/DVN/SF1DRW](https://doi.org/10.7910/DVN/SF1DRW).

## Installation

The R script is written using R 4.0.2. The following R packages are used in this analysis:

### Data handling
- 'tidyverse', including 'readr' (read data), 'tibble' (update dataframes), 'dplyr' (manipulate data), 'ggplot2' (create graphics)
- 'ggalt' (extra coordinate systems, 'geoms', statistical transformations, scales and fonts for 'ggplot2')
- 'gdata' (data manipulation)

### Matching algorithms
- 'MatchIt' (multivariate matching)

### Panel data analysis
- 'plm' (linear models for panel data) 
- 'Hmisc' (miscellaneous functions for data analysis)

### Visualization
- 'ggplot2' (complex plots from data frames)
- 'gridExtra' (miscellaneous functions for "grid" graphics)

## Data Analysis

Original source files have been pre-processed for analysis. To replicate the results, the two datasets are required:

- `PropertyStats.csv`
- `ELC.csv`

The R-code `Albany_teaching_case.R`compares participating properties to a statistical reference group of non-participating properties using statistical techniques. The script implements the propensity score matching procedure to reduce self-selection bias followed by regression adjustments with time and group fixed effects. The choice of matching variables to construct proper counterfactual is at a discretion of students. 
