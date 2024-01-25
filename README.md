
# longterm_synchrony_Cedar_Creek

Code for the paper "Short-term versus multi-decadal responses of
community synchrony, stability, and biodiversity to multiple global
change drivers"

Contains: 1. Raw community data (/data), 2. Cleaned data and subsetted
dataframes (/data_cleaning), 3. Models and analyses (/data_analyses),
and, Supplemental analyses(/data_analyses). See below for additional
detail. Additionally contains functions for calculating the Evar metric
and Box-Cox transformations (/Functions) and exported figure (/Figures) files.

1.  **Data used:**

    *Cedar_Creek_Plant_Taxon_List.csv* - contains data about species
    found at Cedar Creek including fuctional group

    *e001-aboveground-mass-2019-09-13.csv* - contains data from
    long-term nitrogen addition experiment to undisturbed plots

    *e002-aboveground-mass-2019-09-13.csv* -contains data from long-term
    nitrogen addition experiment to disturbed plots

    *final_aboveground_data_for_analyses* - final dataframe used after cleaning the raw
    datasets above and removing certain species for our specific analyses

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

    *spaabundance.csv -* found in /cleaning_CC.R and used to plot species abundance through time

2.  **Data Cleaning**

    *cleaning_CC.R* - script contains data cleaning protocols following
    those found in Seabloom et al. 2020 and DeSiervo et al. 2023 (see
    references for citations).

    *subsetting_CC.R* - script contains data subsetting methods
    intotransient and post-tranisent phases


3.  **Data Analyses**

    *F1_Synchrony_Stability.R* - analyses for global change driver effects on
    synchrony and stability

    *F2_biomass.R - analyses for total biomass through time across
    dominant species*

    *F3_Sync_Stab_Breakdowns.R - analyses for synchrony and stability metric
    breakdown comparisons*

    *F4_S2_Sync_Stab_Relationship.R - analyses for synchrony - stability
    relationship across successional timescales*

    *F5_SEM_Figures_and_Effects.R - analyses for SEM*

    *S1_Control_Treatment_Comparison.R - analyses for comparing two
    control treatments*

    *S3_Comm_Metrics.R - analyses for community metrics across
    successional timescales*

    *S4_SEM_Figures_and_Effects_supp.R - analyses for supplemental SEM*
