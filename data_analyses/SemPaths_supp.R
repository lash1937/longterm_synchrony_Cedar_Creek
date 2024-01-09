############################
###Clean Supplementary SEM code###
############################

#Packages:
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(vegan)

library(lavaan)
library(psych)
library(semPlot)
library(piecewiseSEM)
library(here)

#Read in data and functions
source(here("data_cleaning/subsetting_CC.R"))

###Transformations###----

#Stability
MASS::boxcox(lm(SEM.b.df$Stability ~ 1)) #Determine ideal lambda
SEM.b.df$TStability <- boxcox_transform(SEM.b.df$Stability, -0.3) #Transform var
shapiro.test(SEM.b.df$TStability) #Test for Normality

MASS::boxcox(lm(SEM.a.df$Stability ~ 1))
SEM.a.df$TStability <- boxcox_transform(SEM.a.df$Stability, 0.05)
shapiro.test(SEM.a.df$TStability)

#VR
MASS::boxcox(lm(SEM.b.df$VR ~ 1))
SEM.b.df$TVR <- boxcox_transform(SEM.b.df$VR, 0.5)
shapiro.test(SEM.b.df$TVR)

MASS::boxcox(lm(SEM.a.df$VR ~ 1))
SEM.a.df$TVR <- boxcox_transform(SEM.a.df$VR, 0.4)
shapiro.test(SEM.a.df$TVR)

#Richness
MASS::boxcox(lm(SEM.b.df$Richness ~ 1))
SEM.b.df$TRichness <- boxcox_transform(SEM.b.df$Richness, 0.85)
shapiro.test(SEM.b.df$TRichness)

MASS::boxcox(lm(SEM.a.df$Richness ~ 1))
SEM.a.df$TRichness <- boxcox_transform(SEM.a.df$Richness, -0.05)
shapiro.test(SEM.a.df$TRichness)

#Evenness
MASS::boxcox(lm(SEM.b.df$Evenness ~ 1))
# SEM.b.df$TEvenness <- boxcox_transform(SEM.b.df$Evenness, -0.35)
SEM.b.df$TEvenness <- log(SEM.b.df$Evenness)
shapiro.test(SEM.b.df$TEvenness)

MASS::boxcox(lm(SEM.a.df$Evenness ~ 1))
SEM.a.df$TEvenness <- boxcox_transform(SEM.a.df$Evenness, 0)
shapiro.test(SEM.a.df$TEvenness)

#Mean Biomass
MASS::boxcox(lm(SEM.b.df$mean_biomass ~ 1))
SEM.b.df$Tmean_biomass <- boxcox_transform(SEM.b.df$mean_biomass, 0.2)
shapiro.test(SEM.b.df$Tmean_biomass)
hist(SEM.b.df$Tmean_biomass)

MASS::boxcox(lm(SEM.a.df$mean_biomass ~ 1))
SEM.a.df$Tmean_biomass <- boxcox_transform(SEM.a.df$mean_biomass, 0.1)
shapiro.test(SEM.a.df$Tmean_biomass)
hist(SEM.a.df$Tmean_biomass)

#Community
MASS::boxcox(lm(SEM.b.df$comm ~ 1))
SEM.b.df$Tcomm <- boxcox_transform(SEM.b.df$comm, 0.25)
shapiro.test(SEM.b.df$Tcomm)
hist(SEM.b.df$Tcomm)

MASS::boxcox(lm(SEM.a.df$comm ~ 1))
SEM.a.df$Tcomm <- boxcox_transform(SEM.a.df$comm, 0)
shapiro.test(SEM.a.df$Tcomm)
hist(SEM.a.df$Tcomm)

#Population
MASS::boxcox(lm(SEM.b.df$pop ~ 1))
SEM.b.df$Tpop <- boxcox_transform(SEM.b.df$pop, 0.25)
shapiro.test(SEM.b.df$Tpop)
hist(SEM.b.df$Tpop)

MASS::boxcox(lm(SEM.a.df$pop ~ 1))
SEM.a.df$Tpop <- boxcox_transform(SEM.a.df$pop, -0.15)
shapiro.test(SEM.a.df$Tpop)
hist(SEM.a.df$Tpop)


###SEM path building###----
#Build transient and post-transient models with paths of interest, using transformed variables

#piecewiseSEM model, transient-----
#Add fields in, transient
field.cat <- as.data.frame(model.matrix(~ field, data = SEM.b.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.b.df$uniqueID)

SEM.b.df.cat <- left_join(SEM.b.df, field.cat, by = "uniqueID")
SEM.b.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

#Add fields in, post-transient
field.cat.a <- as.data.frame(model.matrix(~ field, data = SEM.a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.a.df$uniqueID)

SEM.a.df.cat <- left_join(SEM.a.df, field.cat.a, by = "uniqueID")
SEM.a.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

m1 <- 'Tmean_biomass ~ TRichness + TEvenness + Disturbance + Nitrogen + fieldB + fieldC
       Tcomm ~ TRichness + TEvenness + Disturbance + Nitrogen + fieldB + fieldC
       Tpop ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TRichness ~  Nitrogen + Disturbance + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen + Disturbance + fieldB + fieldC
       Tmean_biomass ~~ Tpop
       Tmean_biomass ~~ Tcomm
       Tpop ~~ Tcomm'

#Fit transient years model
m1.fit <- sem(m1, data=SEM.b.df.cat, se="bootstrap", test="bootstrap")
summary(m1.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m1.fit, type="std.all"), 
        file = here::here("data/SEM_supp_transient.rds")) 
  #To read: object <- readRDS(here("data/SEM_supp_transient.rds"))

#Fit post-transient years model
m2.fit <- sem(m1, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
summary(m2.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m2.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m2.fit, type="std.all"), 
        file = here::here("data/SEM_supp_posttransient.rds")) 
#To read: object <- readRDS(here("data/SEM_supp_posttransient.rds"))


# m1psem_supp <- psem(
#   lm(mean_biomass ~ TRichness + TEvenness + Disturbance + Nitrogen + field, data = SEM.b.df),
#   lm(comm ~ TRichness + TEvenness + Disturbance + Nitrogen + field, data = SEM.b.df),
#   lm(pop ~ TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.b.df),
#   lm(TRichness ~  Nitrogen + Disturbance + field, data = SEM.b.df),
#   lm(TEvenness ~  Nitrogen + Disturbance + field, data = SEM.b.df),
#   Tmean_biomass %~~% Tpop,
#   Tmean_biomass %~~% Tcomm,
#   Tpop %~~% Tcomm
# )
# summary(m1psem_supp) 

# bootstrapped model fit
m1psem_boot_supp <- semEff::bootEff(m1psem_supp, R = 4000, parallel = "multicore", ncpus = 4)
m1semeff_supp <- semEff::semEff(m1psem_boot_supp)
summary(m1semeff_supp)
summary(semEff::semEff(m1psem_boot_supp, predictor = "TRichness"), response = "TStability")



#piecewiseSEM model, post-transient------

m2psem_supp <- psem(
  lm(mean_biomass ~ TRichness + TEvenness + Disturbance + Nitrogen + field, data = SEM.a.df),
  lm(comm ~ TRichness + TEvenness + Disturbance + Nitrogen + field, data = SEM.a.df),
  lm(pop ~ TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.a.df),
  lm(TRichness ~  Nitrogen + Disturbance + field, data = SEM.a.df),
  lm(TEvenness ~  Nitrogen + Disturbance + field, data = SEM.a.df),
  mean_biomass %~~% pop,
  mean_biomass %~~% comm,
  pop %~~% comm
)
summary(m2psem_supp) 

# bootstrapped model fit
m2psem_boot_supp <- semEff::bootEff(m2psem_supp, R = 4000, parallel = "multicore", ncpus = 4)
m2semeff_supp <- semEff::semEff(m2psem_boot_supp)
summary(m2semeff)

# save supplemental bootstrapped fits
saveRDS(
  list(
    first7 = m1psem_boot_supp,
    last7 = m2psem_boot_supp
  ),
  file = here("data/bootstrapped_supp_semfits.rds")
)
