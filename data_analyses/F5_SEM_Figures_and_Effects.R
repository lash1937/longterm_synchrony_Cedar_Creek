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
#Calculate indirect effect
ind_eff_NSyS.b <- path.b.coefs[10,3] * path.b.coefs[1,3]
#Create vector of direct paths
NSyS_est.b <- c(
  path.b.coefs[10, 3], #Nitrogen -> Synchrony direct path
  path.b.coefs[1, 3]   #Synchrony -> Stability direct path
)
#Construct covariance matrix
NSyS_cov.b <- matrix(0, 2, 2)
NSyS_cov.b[1, 1] <- path.b.coefs[10, 4]^2         #Synchrony ~ Nitrogen submodel
NSyS_cov.b[2, 2] <- path.b.coefs[1, 4]^2           #Stability ~ Synchrony submodel
#Calculate maximum possible variance
NSyS_cov.b[1, 2] <- path.b.coefs[10, 4] * path.b.coefs[1, 4]
NSyS_cov.b[2,1] <- NSyS_cov.b[1, 2]
#Calculate SE using the delta method
seNSyS.b <- msm::deltamethod(~ x1 * x2, NSyS_est.b, cov = NSyS_cov.b)
#Calculte p-value
pNSyS.b <- pnorm(abs(ind_eff_NSyS.b / seNSyS.b), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
NSyS.b <- c("NSyS", ind_eff_NSyS.b, seNSyS.b, pNSyS.b)

#Disturbance -> Synchrony -> Stability
#Calculate indirect effect
ind_eff_DSyS.b <- path.b.coefs[11,3] * path.b.coefs[1,3]
#Create vector of direct paths
DSyS_est.b <- c(
  path.b.coefs[11, 3], #Dist -> Synchrony direct path
  path.b.coefs[1, 3]   #Synchrony -> Stability direct path
)
#Construct covariance matrix
DSyS_cov.b <- matrix(0, 2, 2)
DSyS_cov.b[1, 1] <- path.b.coefs[11, 4]^2          #Synchrony ~ Dist submodel
DSyS_cov.b[2, 2] <- path.b.coefs[1, 4]^2           #Stability ~ Synchrony submodel
#Calculate maximum possible variance
DSyS_cov.b[1, 2] <- path.b.coefs[11, 4] * path.b.coefs[1, 4]
DSyS_cov.b[2,1] <- DSyS_cov.b[1, 2]
#Calculate SE using the delta method
seDSyS.b <- msm::deltamethod(~ x1 * x2, DSyS_est.b, cov = DSyS_cov.b)
#Calculate p-value
pDSyS.b <- pnorm(abs(ind_eff_DSyS.b / seDSyS.b), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
DSyS.b <- c("DSyS", ind_eff_DSyS.b, seDSyS.b, pDSyS.b)

#Nitrogen -> Richness -> Stability
#Calculate indirect effect
ind_eff_NRS.b <- path.b.coefs[14,3] * path.b.coefs[2,3]
#Create vector of direct paths
NRS_est.b <- c(
  path.b.coefs[14, 3], #Nitrogen -> Richness direct path
  path.b.coefs[2, 3]   #Richness -> Stability direct path
)
#Construct covariance matrix
NRS_cov.b <- matrix(0, 2, 2)
NRS_cov.b[1, 1] <- path.b.coefs[14, 4]^2 #Richness ~ Nitrogen submodel
NRS_cov.b[2, 2] <- path.b.coefs[2, 4]^2  #Stability ~ Richness submodel
#Calculate maximum possible variance
NRS_cov.b[1, 2] <- -path.b.coefs[14, 4] * path.b.coefs[2, 4]
NRS_cov.b[2,1] <- NRS_cov.b[1, 2]
#Calculate SE using the delta method
seNRS.b <- msm::deltamethod(~ x1 * x2, NRS_est.b, cov = NRS_cov.b)
#Calculate p-value
pNRS.b <- pnorm(abs(ind_eff_NRS.b / seNRS.b), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
NRS.b <- c("NRS", ind_eff_NRS.b, seNRS.b, pNRS.b)

#Disturbance -> Richness -> Stability
#Calculate indirect effect
ind_eff_DRS.b <- path.b.coefs[15,3] * path.b.coefs[2,3]
#Create vector of direct paths
DRS_est.b <- c(
  path.b.coefs[15, 3], #Disturbance -> Richness direct path
  path.b.coefs[2, 3]   #Richness -> Stability direct path
)
#Construct covariance matrix
DRS_cov.b <- matrix(0, 2, 2)
DRS_cov.b[1, 1] <- path.b.coefs[15, 4]^2 #Richness ~ Dist submodel
DRS_cov.b[2, 2] <- path.b.coefs[2, 4]^2 #Stability ~ Richness submodel
#Calculate maximum possible variance
DRS_cov.b[1, 2] <- -path.b.coefs[15, 4] * path.b.coefs[2, 4]
DRS_cov.b[2,1] <- DRS_cov.b[1, 2]
#Calculate SE using the delta method
seDRS.b <- msm::deltamethod(~ x1 * x2, DRS_est.b, cov = DRS_cov.b)
#Calculate p-value
pDRS.b <- pnorm(abs(ind_eff_DRS.b / seDRS.b), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
DRS.b <- c("DRS", ind_eff_DRS.b, seDRS.b, pDRS.b)

#Compile results into single dataframe
ind.eff.b <- rbind(NSyS.b, DSyS.b, NRS.b, DRS.b)
colnames(ind.eff.b) <- c("Path", "Effect", "SE", "Pvalue")  



#####Post-transient Phase#####
#Nitrogen -> Synchrony -> Stability
#Calculate indirect effect
ind_eff_NSyS.a <- path.a.coefs[10, 3] * path.a.coefs[1, 3]
#Create vector of direct paths
NSyS_est.a <- c(
  path.a.coefs[10, 3], #Nitrogen -> Synchrony direct path
  path.a.coefs[1, 3]   #Synchrony -> Stability direct path
)
#Construct covariance matrix
NSyS_cov.a <- matrix(0, 2, 2)
NSyS_cov.a[1, 1] <- path.a.coefs[10, 4]^2          #Synchrony ~ Nitrogen submodel
NSyS_cov.a[2, 2] <- path.a.coefs[1, 4]^2           #Stability ~ Synchrony submodel
#Calculate maximum possible variance
NSyS_cov.a[1, 2] <- -path.a.coefs[10, 4] * path.a.coefs[1, 4]
NSyS_cov.a[2,1] <- NSyS_cov.a[1, 2]
#Calculate SE using the delta method
seNSyS.a <- msm::deltamethod(~ x1 * x2, NSyS_est.a, cov = NSyS_cov.a)
#Calculte p-value
pNSyS.a <- pnorm(abs(ind_eff_NSyS.a / seNSyS.a), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
NSyS.a <- c("NSyS", ind_eff_NSyS.a, seNSyS.a, pNSyS.a)

#Disturbance -> Synchrony -> Stability
#Calculate indirect effect
ind_eff_DSyS.a <- path.a.coefs[11, 3] * path.a.coefs[1, 3]
#Create vector of direct paths
DSyS_est.a <- c(
  path.a.coefs[11, 3], #Dist -> Synchrony direct path
  path.a.coefs[1, 3]   #Synchrony -> Stability direct path
)
#Construct covariance matrix
DSyS_cov.a <- matrix(0, 2, 2)
DSyS_cov.a[1, 1] <- path.a.coefs[11, 4]^2 #Synchrony ~ Dist submodel
DSyS_cov.a[2, 2] <- path.a.coefs[1, 4]^2          #Stability ~ Synchrony submodel
#Calculate maximum possible variance
DSyS_cov.a[1, 2] <- path.a.coefs[11, 4] * path.a.coefs[1, 4]
DSyS_cov.a[2,1] <- DSyS_cov.a[1, 2]
#Calculate SE using the delta method
seDSyS.a <- msm::deltamethod(~ x1 * x2, DSyS_est.a, cov = DSyS_cov.a)
#Calculate p-value
pDSyS.a <- pnorm(abs(ind_eff_DSyS.a / seDSyS.a), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
DSyS.a <- c("DSyS", ind_eff_DSyS.a, seDSyS.a, pDSyS.a)

#Nitrogen -> Richness -> Stability
#Calculate indirect effect
ind_eff_NRS.a <- path.a.coefs[14, 3] * path.a.coefs[2, 3]
#Create vector of direct paths
NRS_est.a <- c(
  path.a.coefs[14, 3], #Nitrogen -> Richness direct path
  path.a.coefs[2, 3]   #Richness -> Stability direct path
)
#Construct covariance matrix
NRS_cov.a <- matrix(0, 2, 2)
NRS_cov.a[1, 1] <- path.a.coefs[14, 4]^2 #Richness ~ Nitrogen submodel
NRS_cov.a[2, 2] <- path.a.coefs[2, 4]^2 #Stability ~ Richness submodel
#Calculate maximum possible variance
NRS_cov.a[1, 2] <- -path.a.coefs[14, 4] * path.a.coefs[2, 4]
NRS_cov.a[2,1] <- NRS_cov.a[1, 2]
#Calculate SE using the delta method
seNRS.a <- msm::deltamethod(~ x1 * x2, NRS_est.a, cov = NRS_cov.a)
#Calculate p-value
pNRS.a <- pnorm(abs(ind_eff_NRS.a / seNRS.a), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
NRS.a <- c("NRS", ind_eff_NRS.a, seNRS.a, pNRS.a)

#Disturbance -> Richness -> Stability
#Calculate indirect effect
ind_eff_DRS.a <- path.a.coefs[15,3] * path.a.coefs[2,3]
#Create vector of direct paths
DRS_est.a <- c(
  path.a.coefs[15, 3], #Disturbance -> Richness direct path
  path.a.coefs[2, 3]   #Richness -> Stability direct path
)
#Construct covariance matrix
DRS_cov.a <- matrix(0, 2, 2)
DRS_cov.a[1, 1] <- path.a.coefs[15, 4]^2 #Richness ~ Dist submodel
DRS_cov.a[2, 2] <- path.a.coefs[2, 4]^2 #Stability ~ Richness submodel
#Calculate maximum possible variance
DRS_cov.a[1, 2] <- -path.a.coefs[15, 4] * path.a.coefs[2, 4]
DRS_cov.a[2,1] <- DRS_cov.a[1, 2]
#Calculate SE using the delta method
seDRS.a <- msm::deltamethod(~ x1 * x2, DRS_est.a, cov = DRS_cov.a)
#Calculate p-value
pDRS.a <- pnorm(abs(ind_eff_DRS.a / seDRS.a), lower.tail = F)
#Collect calculated effect, standard error, and p-value together
DRS.a <- c("DRS", ind_eff_DRS.a, seDRS.a, pDRS.a)

#Compile results into single dataframe
ind.eff.a <- matrix(nrow = 4, ncol = 4, 0)
colnames(ind.eff.a) <- c("Path", "Effect", "SE", "Pvalue")  
ind.eff.a[1,] <- NSyS.a
ind.eff.a[2,] <- DSyS.a
ind.eff.a[3,] <- NRS.a
ind.eff.a[4,] <- DRS.a

#Reported indirect effects and errors
ind.eff.b
ind.eff.a

