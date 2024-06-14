#####################################
# This script produces summary outputs of structural equation 
# models used for creating supplemental SEM pathway figures. These figures 
# are equivlent to the main text Figure 5, but looks at relationship pathways
# when all years are included in the SEM.
#####################################

#Packages:
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(ggpubr)

library(lavaan)
library(msm)
library(psych)
library(here)

#Read in data and functions
source(here::here("data_cleaning/subsetting_CC.R"))

###Data Transformations###----
#Nitrogen
SEM.df$logN <- log(SEM.df$Nitrogen + 1)

#Stability
MASS::boxcox(lm(SEM.df$Stability ~ 1)) #Determine ideal lambda
SEM.df$TStability <- boxcox_transform(SEM.df$Stability, -0.2) #Transform variable
shapiro.test(SEM.df$TStability) #Test for Normality

#VR
MASS::boxcox(lm(SEM.df$VR ~ 1))
SEM.df$TVR <- boxcox_transform(SEM.df$VR, 0.3)
shapiro.test(SEM.df$TVR)

#Richness
MASS::boxcox(lm(SEM.df$Richness ~ 1))
SEM.df$TRichness <- boxcox_transform(SEM.df$Richness, -0.05)
shapiro.test(SEM.df$TRichness)

#Evenness
MASS::boxcox(lm(SEM.df$Evenness ~ 1))
SEM.df$TEvenness <- boxcox_transform(SEM.df$Evenness, -0.4)
shapiro.test(SEM.df$TEvenness)


############################
###SEM model building###
# Produce summary outputs of direct effects of exogenous factors
# on community properties across all years in the dataset
############################

#Build model used in SEM
#Set insignificant pathways to 0, let significant pathways freely vary
m1supp <- 'TStability ~ TVR + c(0, 0)*TRichness + c(NA, 0)*TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ c(0, NA)*TRichness + c(0, 0)*TEvenness + c(0, NA)*Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ c(NA, 0)*TRichness + c(NA, 0)*Nitrogen  + fieldB + fieldC'

#Fit ALL years model 
msupp.fit <- sem(m1supp, data=SEM.df, group = "Disturbance")
summary(msupp.fit, fit.measures=TRUE, stand=TRUE, rsq=TRUE)
standardizedSolution(msupp.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(msupp.fit, type="std.all"),
        file = here::here("data/SEM_allyears.rds"))
#object <- readRDS(here("data/SEM_allyears.rds"))