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
library(ggpubr)
library(vegan)

library(piecewiseSEM)
library(nlme)
library(msm)
library(psych)
library(here)

#Read in data and functions
source(here::here("data_cleaning/subsetting_CC.R"))

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

MASS::boxcox(lm(SEM.a.df$mean_biomass ~ 1))
SEM.a.df$Tmean_biomass <- boxcox_transform(SEM.a.df$mean_biomass, 0.1)
shapiro.test(SEM.a.df$Tmean_biomass)

#Community Variance
MASS::boxcox(lm(SEM.b.df$comm ~ 1))
SEM.b.df$Tcomm <- boxcox_transform(SEM.b.df$comm, 0.25)
shapiro.test(SEM.b.df$Tcomm)

MASS::boxcox(lm(SEM.a.df$comm ~ 1))
SEM.a.df$Tcomm <- boxcox_transform(SEM.a.df$comm, 0)
shapiro.test(SEM.a.df$Tcomm)

#Population Variance
MASS::boxcox(lm(SEM.b.df$pop ~ 1))
SEM.b.df$Tpop <- boxcox_transform(SEM.b.df$pop, 0.25)
shapiro.test(SEM.b.df$Tpop)

MASS::boxcox(lm(SEM.a.df$pop ~ 1))
SEM.a.df$Tpop <- boxcox_transform(SEM.a.df$pop, -0.15)
shapiro.test(SEM.a.df$Tpop)


##Scale all non-categorical variables
SEM.b.df <- SEM.b.df %>% 
  mutate_at(c('Nitrogen', 'TStability', 'TVR', 'TRichness', 'TEvenness',
              'Tmean_biomass', 'Tcomm', 'Tpop'), ~(scale(.)
                                                  %>% as.vector))

SEM.a.df <- SEM.a.df %>% 
  mutate_at(c('Nitrogen', 'TStability', 'TVR', 'TRichness', 'TEvenness',
              'Tmean_biomass', 'Tcomm', 'Tpop'), ~(scale(.)
                                                  %>% as.vector))


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

###Piecewise model, transient phase
#Model structure
m1psem <- psem(
  nlme::lme(Tmean_biomass ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + 
              fieldC, random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(Tcomm ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + fieldC,
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(Tpop ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + fieldC, 
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(TRichness ~ Nitrogen + Disturbance + fieldB + fieldC, 
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(TEvenness ~  Nitrogen + Disturbance + TRichness + fieldB + fieldC,
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  data = SEM.b.df
)
#Add in partial correlations between decomposed synchrony and stability components 
m1psem <- update(m1psem, Tmean_biomass %~~% Tpop, Tmean_biomass %~~% Tcomm,
                           Tcomm %~~% Tpop)

#Report path coefficients, standard errors, and p-values
#for transient phase data
coefs(m1psem)


###Piecewise model, post-transient phase
#Model structure
m2psem <- psem(
  nlme::lme(Tmean_biomass ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + 
              fieldC, random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(Tcomm ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + fieldC,
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(Tpop ~ Nitrogen + Disturbance + TRichness + TEvenness + fieldB + fieldC, 
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(TRichness ~ Nitrogen + Disturbance + fieldB + fieldC, 
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(TEvenness ~  Nitrogen + Disturbance + TRichness + fieldB + fieldC,
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  data = SEM.a.df
)
#Add in partial correlations between decomposed synchrony and stability components 
m2psem <- update(m2psem, Tmean_biomass %~~% Tpop, Tmean_biomass %~~% Tcomm,
                 Tcomm %~~% Tpop)

#Report path coefficients, standard errors, and p-values
#for post-transient phase data
coefs(m2psem)
