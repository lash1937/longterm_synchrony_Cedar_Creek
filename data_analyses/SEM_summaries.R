################################################################
# With this script, we record code used to summarize the results
# from the SEMs
################################################################

##### Document setup #####
  library(semEff)
  library(here)
  source(here::here("R/functions.R"))

##### Load main figure data #####
  sems_boot <- readRDS(here("data/bootstrapped_semfits.rds"))

##### Summaries #####
  semeff_f7_full <- semEff(sems_boot$first7)
  semeff_l7_full <- semEff(sems_boot$last7)

  #Without spp files (using new data - da.widesynch)
  export_sem_table(semeff_f7_full, outfile = here("Tables/first7_sem_results.tex"), longtable = T)
  export_sem_table(semeff_l7_full, outfile = here("Tables/last7_sem_results.tex"), longtable = T)
  #With spp files (using old data - da.widewoodysynch)
  # export_sem_table(semeff_f7_full, outfile = here("Tables/woody_first7_sem_results.tex"), longtable = T)
  # export_sem_table(semeff_l7_full, outfile = here("Tables/woody_last7_sem_results.tex"), longtable = T)

  summary(semEff(sems_boot$first7, predictors = "TRichness"), response = "TStability")  
  
##### Load supplemental data #####
  sems_boot_supp <- readRDS(here("data/bootstrapped_supp_semfits.rds"))
  
##### Summaries #####
  semeff_f7_full_supp <- semEff(sems_boot_supp$first7)
  semeff_l7_full_supp <- semEff(sems_boot_supp$last7)
  
  export_sem_table(semeff_f7_full_supp, outfile = here("Tables/supp_first7_sem_results.tex"), longtable = T)
  export_sem_table(semeff_l7_full_supp, outfile = here("Tables/supp_last7_sem_results.tex"), longtable = T)
  
  summary(semEff(sems_boot_supp$first7, predictors = "TRichness"), response = "TStability")  
  