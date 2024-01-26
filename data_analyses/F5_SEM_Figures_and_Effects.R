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

library(piecewiseSEM)
library(nlme)
library(msm)
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



##Scale all non-categorical variables
SEM.b.df <- SEM.b.df %>% 
  mutate_at(c('Nitrogen', 'TStability', 'TVR', 'TRichness', 'TEvenness'), ~(scale(.)
                                                                            %>% as.vector))

SEM.a.df <- SEM.a.df %>% 
  mutate_at(c('Nitrogen', 'TStability', 'TVR', 'TRichness', 'TEvenness'), ~(scale(.)
                                                                            %>% as.vector))

############################
###SEM model building###
# Produce summary outputs of direct effects of exogenous factors
# on community properties
############################

###Piecewise model, transient phase
#Model structure
m1psem <- psem(
  nlme::lme(TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + fieldB + 
              fieldC, random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC,
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  nlme::lme(TRichness ~  Nitrogen + Disturbance + fieldB + fieldC, random = (~1|grid), 
            data = SEM.b.df, method = "ML"),
  nlme::lme(TEvenness ~  TRichness + Nitrogen + Disturbance + fieldB + fieldC,
            random = (~1|grid), data = SEM.b.df, method = "ML"),
  data = SEM.b.df
)

#Report path coefficients, standard errors, and p-values
#for transient phase data
coefs(m1psem)

#Calculate Confidence Intervals, sub-models 1-4
intervals(m1psem[[1]], which = "fixed") #Effects of variables on Stability
intervals(m1psem[[2]], which = "fixed") #Effects of variables on Synchrony
intervals(m1psem[[3]], which = "fixed") #Effects of variables on Richness
intervals(m1psem[[4]], which = "fixed") #Effects of variables on Richness


###Piecewise model, post-transient phase
#Model structure
m2psem <- psem(
  nlme::lme(TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + fieldB + 
              fieldC, random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC,
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  nlme::lme(TRichness ~  Nitrogen + Disturbance + fieldB + fieldC, random = (~1|grid), 
            data = SEM.a.df, method = "ML"),
  nlme::lme(TEvenness ~  TRichness + Nitrogen + Disturbance + fieldB + fieldC,
            random = (~1|grid), data = SEM.a.df, method = "ML"),
  data = SEM.a.df
)

#Report path coefficients, standard errors, and p-values
#for post-transient phase data
coefs(m2psem)

#Calculate Confidence Intervals, sub-models 1-4
intervals(m2psem[[1]], which = "fixed") #Effects of variables on Stability
intervals(m2psem[[2]], which = "fixed") #Effects of variables on Synchrony
intervals(m2psem[[3]], which = "fixed") #Effects of variables on Richness
intervals(m2psem[[4]], which = "fixed") #Effects of variables on Richness


############################
###(In)Direct Effects and Mediation###
# Produce summary outputs of selected indirect effects 
# of exogenous factors on community properties
############################

#Save path coefficients for calculation of indirect effects
path.b.coefs <- coefs(m1psem)
path.a.coefs <- coefs(m2psem)

#####Transient Phase#####
#Nitrogen -> Synchrony -> Stability
ind_eff_NSyS <- coefs(m1psem)[10,3] * coefs(m1psem)[1,3] #Indirect effect

NSyS_est <- c(
  path.b.coefs[10, 3], #Nitrogen -> Synchrony direct path
  path.b.coefs[1, 3]   #Synchrony -> Stability direct path
)
#Construct covariance matrix
NSyS_cov <- matrix(0, 2, 2)
NSyS_cov[1, 1] <- vcov(m1psem[[2]])["Nitrogen", "Nitrogen"] #Synchrony ~ Nitrogen submodel
NSyS_cov[2, 2] <- vcov(m1psem[[1]])["TVR", "TVR"]           #Stability ~ Synchrony submodel

#Calculate maximum possible variance
NSyS_cov[1, 2] <- sqrt(NSyS_cov[1,1]) * sqrt(NSyS_cov[2,2])
NSyS_cov[2,1] <- NSyS_cov[1, 2]

#Calculate SE using the delta method
seNSyS <- msm::deltamethod(~ x1 * x2, NSyS_est, cov = NSyS_cov)

# the z-stat is just the estimate divided by its SE
# when the null hypothesis is that the effect is zero
# the z-stat can be compared to a standard normal to get the p-value
# for the test
pNSyS <- pnorm(abs(ind_eff_NSyS / seNSyS), lower.tail = F)

#Collect calculated effect, standard error, and p-value together
NSyS <- c("NSyS", ind_eff_NSyS, seNSyS, pNSyS)

#Disturbance -> Synchrony -> Stability

#Nitrogen -> Richness -> Stability

#Disturbance -> Richness -> Stability

#####Post-transient Phase#####
#Nitrogen -> Synchrony -> Stability

#Disturbance -> Synchrony -> Stability

#Nitrogen -> Richness -> Stability

#Disturbance -> Richness -> Stability

#Compile into single dataframe
ind.eff.m <- matrix(nrow = 4, ncol = 4, 0)
colnames(ind.eff.m) <- c("Path", "Effect", "SE", "Pvalue")  
ind.eff.m[1,] <- NSyS
ind.eff.m[2,] <- DSyS
ind.eff.m[3,] <- NRS
ind.eff.m[4,] <- DRS

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
'

