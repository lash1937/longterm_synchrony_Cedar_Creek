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

#Create model structure
#Build model used in all SEMs
m1 <- 'TStability ~ TVR + TRichness + TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ TRichness + TEvenness + Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen  + fieldB + fieldC'

#Fit ALL years model 
m3.fit <- sem(m1, data=SEM.df, group = "Disturbance")
#add in bootstrapping when model is determined to be fine se="bootstrap", test="bootstrap"
summary(m3.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m3.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(m3.fit, type="std.all"),
        file = here::here("data/SEM_allyears.rds"))
object <- readRDS(here("data/SEM_allyears.rds"))
