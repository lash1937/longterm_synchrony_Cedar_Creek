#####################################
# This script produces summary outputs of structural equation 
# models used for creating supplemental SEM pathway figures. These figures 
# break down the variance ratio and inverse coefficient of variation into
# their mathematical components, and report all direct effects.
#####################################

#Packages:
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(vegan)

library(lavaan)
library(psych)
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

#Community Variance
MASS::boxcox(lm(SEM.b.df$comm ~ 1))
SEM.b.df$Tcomm <- boxcox_transform(SEM.b.df$comm, 0.25)
shapiro.test(SEM.b.df$Tcomm)
hist(SEM.b.df$Tcomm)

MASS::boxcox(lm(SEM.a.df$comm ~ 1))
SEM.a.df$Tcomm <- boxcox_transform(SEM.a.df$comm, 0)
shapiro.test(SEM.a.df$Tcomm)
hist(SEM.a.df$Tcomm)

#Population Variance
MASS::boxcox(lm(SEM.b.df$pop ~ 1))
SEM.b.df$Tpop <- boxcox_transform(SEM.b.df$pop, 0.25)
shapiro.test(SEM.b.df$Tpop)
hist(SEM.b.df$Tpop)

MASS::boxcox(lm(SEM.a.df$pop ~ 1))
SEM.a.df$Tpop <- boxcox_transform(SEM.a.df$pop, -0.15)
shapiro.test(SEM.a.df$Tpop)
hist(SEM.a.df$Tpop)


############################
###SEM model building###
# Produce summary outputs of direct effects of exogenous factors
# on community properties
############################

#Create model structure
m1 <- 'Tmean_biomass ~ TRichness + TEvenness + Disturbance + Nitrogen + fieldB + fieldC
       Tcomm ~ TRichness + TEvenness + Disturbance + Nitrogen + fieldB + fieldC
       Tpop ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TRichness ~  Nitrogen + Disturbance + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen + Disturbance + fieldB + fieldC
       Tmean_biomass ~~ Tpop
       Tmean_biomass ~~ Tcomm
       Tpop ~~ Tcomm'

###Lavaan model, transient phase
#Add in field to dataframe
field.cat <- as.data.frame(model.matrix(~ field, data = SEM.b.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.b.df$uniqueID)
SEM.b.df.cat <- left_join(SEM.b.df, field.cat, by = "uniqueID")
SEM.b.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

#Fit model to transient phase data
m1.fit <- sem(m1, data=SEM.b.df.cat, se="bootstrap", test="bootstrap")
standardizedSolution(m1.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m1.fit, type="std.all"), 
        file = here::here("data/SEM_supp_transient.rds")) 

#Add in field to dataframe
field.cat.a <- as.data.frame(model.matrix(~ field, data = SEM.a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.a.df$uniqueID)
SEM.a.df.cat <- left_join(SEM.a.df, field.cat.a, by = "uniqueID")
SEM.a.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

#Fit model to post-transient phase data
m2.fit <- sem(m1, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
standardizedSolution(m2.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m2.fit, type="std.all"), 
        file = here::here("data/SEM_supp_posttransient.rds"))