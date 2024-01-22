---
editor_options: 
  markdown: 
    wrap: 72
---

# longterm_synchrony_Cedar_Creek

Code for the paper "Short-term versus multi-decadal responses of
community synchrony, stability, and biodiversity to multiple global
change drivers"

Contains: 1. Raw community data (/data), 2. Cleaned data and subsetted
dataframes (/data_cleaning), 3. Models and analyses (/data_analyses),
and, Supplemental analyses(/data_analyses). See below for additional
detail. Additionally contains functions for calculating the Evar metric
and Box-Cox transformations (/Functions) and exported figure (/Figures)
and table (/Tables) files.

1.  **Data used:**

    *Cedar_Creek_Plant_Taxon_List.csv* - contains data about species
    found at Cedar Creek including fuctional group

    *e001-aboveground-mass-2019-09-13.csv* - contains data from
    long-term nitrogen addition experiment to undisturbed plots

    *e002-aboveground-mass-2019-09-13.csv* -contains data from long-term
    nitrogen addition experiment to disturbed plots

    *NEW.sp.decisions* - final dataframe used after cleaning the raw
    datasets above and removing species for our specific analyses

    **Data listed but not used:**

    *e001-belowground-mass-2019-09-13.csv -* found in /cleaning_CC.R
    where data cleaning protocols follow those found in Seabloom et al.
    2020 and DeSiervo et al. 2023 (see references for citations).

    *E001-E002-Soil-CN-2018 -* ound in /cleaning_CC.R where data
    cleaning protocols follow those found in Seabloom et al. 2020 and
    DeSiervo et al. 2023 (see references for citations).

    *e001-soil-cn-2019-09-13.csv - f*ound in /cleaning_CC.R where data
    cleaning protocols follow those found in Seabloom et al. 2020 and
    DeSiervo et al. 2023 (see references for citations).

    *sp.decisions.csv -* found in /cleaning_CC.R used to test
    sensitivity of removing woody species

    *spaabundance.csv -* found in /cleaning_CC.R and used to plot

2.  **Data Cleaning**

    *cleaning_CC.R* - script contains data cleaning protocols following
    those found in Seabloom et al. 2020 and DeSiervo et al. 2023 (see
    references for citations).

    *subsetting_CC.R* - script contains data subsetting methods
    intotransient and post-tranisent phases

<!-- -->

3.  **Data Analyses**

    *Clean_figure_1.R* - analyses for global change driver effects on
    synchrony and stability

    *Clean_figure_2.R - analyses for total biomass through time across
    dominant species*

    *Clean_figure_3.R - analyses for synchrony and stability metric
    breakdown comparisons*

    *Clean_figure_4\_S2.R - analyses for synchrony - stability
    relationship across successional timescales*

    *5_SEM_Figures_and_Effects.R - analyses for SEM*

    *Control_treatment_comparison_figS1.R - analyses for comparing two
    control treatments*

    *Clean_supp_figS3.R - analyses for community metrics across
    successional timescales*

    *S4_SEM_Figures_and_Effects_supp.R - analyses for supplemental SEM*
