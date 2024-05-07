#####################################
# This script produces summary outputs of structural equation models,
# which are used for creating SEM pathway figures. These figures 
# represent the effects of nitrogen addition, dependant on the 
# presence of soil disturbance treatments, on community properties, 
# and report all direct and selected indirect effects.
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
SEM.b.df$logN <- log(SEM.b.df$Nitrogen + 1)
SEM.a.df$logN <- log(SEM.a.df$Nitrogen + 1)

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

#Build model used in transient SEMs
m1 <- 'TStability ~ TVR + TRichness + c(0,NA)*TEvenness + c(NA,0)*Nitrogen + fieldB + fieldC
       TVR ~ c(0, NA)*TRichness + c(0,0)*TEvenness + Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ c(0,0)*TRichness + c(0, NA)*Nitrogen  + fieldB + fieldC'



###Lavaan model, transient phase
m1.fit <- sem(m1, data=SEM.b.df, group = "Disturbance")
summary(m1.fit, fit.measures=TRUE, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(m1.fit, type="std.all"),
        file = here::here("data/SEM_transient.rds"))
#object <- readRDS(here("data/SEM_transient.rds"))



###Lavaan model, post-transient phase
m2.fit <- sem(m1, data=SEM.a.df, group = "Disturbance", se="bootstrap", test="bootstrap")
summary(m2.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m2.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(m2.fit, type="std.all"),
        file = here::here("data/SEM_posttransient.rds"))
#object <- readRDS(here("data/SEM_posttransient.rds"))


#########################################################################################
######Testing unsaturated model#######

#c(nodist, dist) order of setting paths to be zero

#Unsaturated model tests - play with this:
m.sat.trans <- 'TStability ~ TVR + TRichness + c(0,NA)*TEvenness + c(NA,0)*Nitrogen + fieldB + fieldC
       TVR ~ c(0, NA)*TRichness + c(0,0)*TEvenness + Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ c(0,0)*TRichness + c(0, NA)*Nitrogen  + fieldB + fieldC'

m.sat.posttrans <- 'TStability ~ TVR + c(0,0)*TRichness + c(0,0)*TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ c(0,0)*TRichness + c(0,0)*TEvenness + c(NA, 0)*Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen  + fieldB + fieldC'

#Transient regression equal
m.sat.regequal <- 'TStability ~ TVR + TRichness + TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ c(0,0)*TRichness + c(0,0)*TEvenness + Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ c(0,0)*TRichness + Nitrogen  + fieldB + fieldC'

#PostTransient regression equal
m.sat.regequal <- 'TStability ~ TVR + TRichness + TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ c(0,0)*TRichness + c(0,0)*TEvenness + Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ c(0,0)*TRichness + Nitrogen  + fieldB + fieldC'

m.sat.posttrans.regequal <- 'TStability ~ TVR + c(0,0)*TRichness + c(0,0)*TEvenness + Nitrogen + fieldB + fieldC
       TVR ~ c(0,0)*TRichness + c(0,0)*TEvenness + c(0,0)*Nitrogen  + fieldB + fieldC
       TRichness ~  Nitrogen  + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen  + fieldB + fieldC'

###Lavaan model, transient phase
m3.fit <- sem(m.sat.trans, data=SEM.b.df, group = "Disturbance")
summary(m3.fit, fit.measures=TRUE, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")

###Lavaan model, transient phase
m4.fit <- sem(m.sat.posttrans, data=SEM.a.df, group = "Disturbance")
summary(m4.fit, fit.measures=TRUE, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")

###Lavaan model, transient phase - equal group regression
m3.fit.equalregression <- sem(m.sat.regequal, data=SEM.b.df, group = "Disturbance", 
              group.equal = "regressions")
summary(m3.fit, fit.measures=TRUE, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")

lavTestLRT(m3.fit, m3.fit.equalregression)

###Lavaan model, posttransient phase - equal group regression
m4.fit.equalregression <- sem(m.sat.posttrans.regequal, data=SEM.a.df, group = "Disturbance", 
                              group.equal = "regressions")
summary(m4.fit.equalregression, fit.measures=TRUE, stand=TRUE, rsq=TRUE)

lavTestLRT(m4.fit, m4.fit.equalregression)









############################
###(In)Direct Effects and Mediation###
# Produce summary outputs of selected indirect effects 
# of exogenous factors on community properties
############################

#Set vectors, so that each model (2 total) has it's own independent indirect effect#####

m.indirect <- '#Direct effects on Stab
                TStability ~ c(c1.g1, c1.g2)*Nitrogen + c(b2.g1, b2.g2)*TEvenness + 
                c(b4.g1, b4.g2)*TRichness + c(d1.g1, d1.g2)*TVR
               #Mediators
                TVR ~ c(a2.g1, a2.g2)*Nitrogen + c(b3.g1, b3.g2)*TRichness + 
                c(b1.g1, b1.g2)*TEvenness
                TRichness ~ c(a3.g1, a3.g2)*Nitrogen  
                TEvenness ~ c(e1.g1, e1.g2)*TRichness + c(a1.g1, a1.g2)*Nitrogen 
               #Extra variables
                TStability ~ fieldB + fieldC
                TVR ~ fieldB + fieldC
                TRichness ~ fieldB + fieldC
                TEvenness ~ fieldB + fieldC
                
               #Indirect effects Group 1
                g1ind_eff_NSyS := a2.g1 * d1.g1
                g1ind_eff_NRS := a3.g1 * b4.g1
                
                #Indirect effects Group 2
                g2ind_eff_NSyS := a2.g2 * d1.g2
                g2ind_eff_NRS := a3.g2 * b4.g2
'

#Fit before years model
m.indirect.fit.b <- sem(m.indirect, data=SEM.b.df, group = "Disturbance")
summary(m.indirect.fit.b, stand=TRUE, rsq=TRUE) #look at rsq values
standardizedSolution(m.indirect.fit.b, type="std.all")

#Save data
saveRDS(standardizedSolution(m.indirect.fit.b, type="std.all"), 
        file = here::here("data/SEM_indirect_transient.rds")) 

#Fit after years model
m.indirect.fit.a <- sem(m.indirect, data=SEM.a.df, se="bootstrap", test="bootstrap", 
                        se="bootstrap", test="bootstrap")
summary(m.indirect.fit.a, stand=TRUE, rsq=TRUE)
standardizedSolution(m.indirect.fit.a, type="std.all")

#Save data
saveRDS(standardizedSolution(m.indirect.fit.a, type="std.all"), 
        file = here::here("data/SEM_indirect_posttransient.rds")) 



######################
# Compare SEM pathway outputs when using a 10 year window
# instead of a 7 year window

###Data Transformations###----
#Nitrogen
SEM.10b.df$logN <- log(SEM.10b.df$Nitrogen + 1)
SEM.10a.df$logN <- log(SEM.10a.df$Nitrogen + 1)

#Stability
MASS::boxcox(lm(SEM.10b.df$Stability ~ 1)) #Determine ideal lambda
SEM.10b.df$TStability <- boxcox_transform(SEM.10b.df$Stability, -0.4) #Transform variable
shapiro.test(SEM.10b.df$TStability) #Test for Normality

MASS::boxcox(lm(SEM.10a.df$Stability ~ 1))
SEM.10a.df$TStability <- boxcox_transform(SEM.10a.df$Stability, 0.05)
shapiro.test(SEM.10a.df$TStability)

#VR
MASS::boxcox(lm(SEM.10b.df$VR ~ 1))
SEM.10b.df$TVR <- boxcox_transform(SEM.10b.df$VR, 0.6)
shapiro.test(SEM.10b.df$TVR)

MASS::boxcox(lm(SEM.10a.df$VR ~ 1))
SEM.10a.df$TVR <- boxcox_transform(SEM.10a.df$VR, 0.5)
shapiro.test(SEM.10a.df$TVR)

#Richness
MASS::boxcox(lm(SEM.10b.df$Richness ~ 1))
SEM.10b.df$TRichness <- boxcox_transform(SEM.10b.df$Richness, 0.5)
shapiro.test(SEM.10b.df$TRichness)

MASS::boxcox(lm(SEM.10a.df$Richness ~ 1))
SEM.10a.df$TRichness <- boxcox_transform(SEM.10a.df$Richness, -0.02)
shapiro.test(SEM.10a.df$TRichness)

#Evenness
MASS::boxcox(lm(SEM.10b.df$Evenness ~ 1))
SEM.10b.df$TEvenness <- boxcox_transform(SEM.10a.df$Evenness, -0.2)
shapiro.test(SEM.10b.df$TEvenness)

MASS::boxcox(lm(SEM.10a.df$Evenness ~ 1))
SEM.10a.df$TEvenness <- boxcox_transform(SEM.10a.df$Evenness, -0.3)
shapiro.test(SEM.10a.df$TEvenness)


###Lavaan model, 10 year transient phase
m4.fit <- sem(m1, data=SEM.10b.df, group = "Disturbance", se="bootstrap", test="bootstrap")
summary(m4.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m4.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(m4.fit, type="std.all"),
        file = here::here("data/SEM_10yr_transient.rds"))
#object <- readRDS(here("data/SEM_10yr_transient.rds"))

###Lavaan model, 10 year post-transient phase
m5.fit <- sem(m1, data=SEM.10a.df, group = "Disturbance", se="bootstrap", test="bootstrap")
summary(m5.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m5.fit, type="std.all")

#Save data
saveRDS(standardizedSolution(m5.fit, type="std.all"),
        file = here::here("data/SEM_10yr_posttransient.rds"))
#object <- readRDS(here("data/SEM_10yr_posttransient.rds"))

