#####################################
# This script produces summary outputs of structural equation models,
# which are used for creating SEM pathway figures. These figures 
# represent the effects of nitrogen addition and soil 
# disturbance on community properties, and report all 
# direct and selected indirect effects.
#####################################

#Packages:
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(ggpubr)
library(vegan)

library(lavaan)
library(psych)
library(here)

#Read in data and functions
source(here::here("data_cleaning/subsetting_CC.R"))

###Data Transformations###----
#Stability
MASS::boxcox(lm(SEM.b.df$Stability ~ 1)) #Determine ideal lambda
SEM.b.df$TStability <- boxcox_transform(SEM.b.df$Stability, -0.3) #Transform variable
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
SEM.b.df$TEvenness <- log(SEM.b.df$Evenness)

MASS::boxcox(lm(SEM.a.df$Evenness ~ 1))
SEM.a.df$TEvenness <- boxcox_transform(SEM.a.df$Evenness, 0)
shapiro.test(SEM.a.df$TEvenness)

############################
###SEM model building###
# Produce summary outputs of direct effects of exogenous factors
# on community properties
############################

#Create model structure
m1 <- 'TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TRichness ~  Nitrogen + Disturbance + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen + Disturbance + fieldB + fieldC'

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
#Save data summary output
saveRDS(standardizedSolution(m1.fit, type="std.all"), 
        file = here::here("data/SEM_transient.rds")) 

###Lavaan model, post-transient phase
#Add in field to dataframe
field.cat.a <- as.data.frame(model.matrix(~ field, data = SEM.a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.a.df$uniqueID)
SEM.a.df.cat <- left_join(SEM.a.df, field.cat.a, by = "uniqueID")
SEM.a.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

#Fit model to post-transient phase data
m2.fit <- sem(m1, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
standardizedSolution(m2.fit, type="std.all")
#Save data summary output
saveRDS(standardizedSolution(m2.fit, type="std.all"), 
        file = here::here("data/SEM_posttransient.rds")) 

############################
###(In)Direct Effects and Mediation###
# Produce summary outputs of selected indirect effects 
# of exogenous factors on community properties
############################

m.indirect <- '#Direct effects on Stability
                TStability ~ c1*Nitrogen + c2*Disturbance+ b2*TEvenness + b4*TRichness + 
                              d1*TVR
               #Direct effects on Mediators
                TVR ~ a3*Nitrogen + a6*Disturbance + b3*TRichness + b1*TEvenness
                TRichness ~ a2*Nitrogen + a5*Disturbance
                TEvenness ~ e1*TRichness + a1*Nitrogen + a4*Disturbance
               #Extra variables
                TStability ~ fieldB + fieldC
                TVR ~ fieldB + fieldC
                TRichness ~ fieldB + fieldC
                TEvenness ~ fieldB + fieldC
               #Indirect effects
                ind_eff_NSyS := a3 * d1
                ind_eff_DSyS := a6 * d1
                ind_eff_NRS := a2 * b4
                ind_eff_DRS := a5 * b4
               #Total effects
                total_eff_N := c1 + a1*b1*d1 + a1*b2+ a2*e1*b1*d1 + a2*e1*b2 + a2*b3*d1 + 
                           a2*b4 + a3*d1 
                total_eff_D := c2 + a5*b3*d1 + a5*b4 + a4*b1*d1 + a5*e1*b1*d1 + a5*b3*d1 +
                           a4*b2 + a6*d1
                total_eff := c1 + a1*b1*d1 + a1*b2+ a2*e1*b1*d1 + a2*e1*b2 + a2*b3*d1 + 
                           a2*b4 + a3*d1 + c2 + a5*b3*d1 + a5*b4 + a4*b1*d1 + a5*e1*b1*d1 + 
                           a5*b3*d1 + a4*b2 + a6*d1
'

#Fit model to transient phase data
m.indirect.fit.b <- sem(m.indirect, data=SEM.b.df.cat, se="bootstrap", test="bootstrap")
standardizedSolution(m.indirect.fit.b, type="std.all")
#Save data
saveRDS(standardizedSolution(m.indirect.fit.b, type="std.all"), 
        file = here::here("data/SEM_indirect_transient.rds")) 

#Fit model to post-transient phase data
m.indirect.fit.a <- sem(m.indirect, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
standardizedSolution(m.indirect.fit.a, type="std.all")
#Save data
saveRDS(standardizedSolution(m.indirect.fit.a, type="std.all"), 
        file = here::here("data/SEM_indirect_posttransient.rds")) 
