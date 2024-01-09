############################
###Clean SEM code###
############################

#Packages:
library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(ggpubr)
library(vegan)
library(nlme)

library(lavaan)
library(psych)
library(semPlot)
library(piecewiseSEM)
library(lme4)
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
SEM.b.df$TEvenness <- log(SEM.b.df$Evenness)

MASS::boxcox(lm(SEM.a.df$Evenness ~ 1))
SEM.a.df$TEvenness <- boxcox_transform(SEM.a.df$Evenness, 0)
shapiro.test(SEM.a.df$TEvenness)

###SEM path building###----
#Build transient and post-transient models with paths of interest, using transformed variables

###Lavaan model
#Add in field as separate columns
field.cat <- as.data.frame(model.matrix(~ field, data = SEM.b.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.b.df$uniqueID)

SEM.b.df.cat <- left_join(SEM.b.df, field.cat, by = "uniqueID")
SEM.b.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

m1 <- 'TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + fieldB + fieldC
       TRichness ~  Nitrogen + Disturbance + fieldB + fieldC
       TEvenness ~ TRichness + Nitrogen + Disturbance + fieldB + fieldC'

#Fit before years model
m1.fit <- sem(m1, data=SEM.b.df.cat, se="bootstrap", test="bootstrap")
summary(m1.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m1.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m1.fit, type="std.all"), 
        file = here::here("data/SEM_transient.rds")) 
object <- readRDS(here("data/SEM_transient.rds"))

#Code for visualizing rudimentary pathways
# semPaths(m1.fit, what="std", nCharNodes = 0, sizeMan = 8, color = "lightblue", layout = "tree",
#          residuals = FALSE, edge.label.cex = 1.5, label.font = 2, label.cex = 2, exoCov = FALSE,
#          curvePivot = TRUE, sizeLat2 = 5, fade = FALSE, edge.label.position = .35)

# #piecewiseSEM model, transient-----
# m1psem <- psem(
#   lm(TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.b.df),
#   lm(TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.b.df),
#   lm(TRichness ~  Nitrogen + Disturbance + field, data = SEM.b.df),
#   lm(TEvenness ~  TRichness + Nitrogen + Disturbance + field, data = SEM.b.df)
# )
# #Path coefficient and significance summary
# summary(m1psem)
# 
# #bootstrapped model fit for transient model
# m1psem_boot <- semEff::bootEff(m1psem, R = 4000, parallel = "multicore", ncpus = 4)
# m1semeff <- semEff::semEff(m1psem_boot)
# #Effects table summary and significance
# summary(m1semeff)
# summary(semEff::semEff(m1psem_boot, predictor = "TRichness"), response = "TStability")



#piecewiseSEM model, post-transient------
field.cat.a <- as.data.frame(model.matrix(~ field, data = SEM.a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.a.df$uniqueID)

SEM.a.df.cat <- left_join(SEM.a.df, field.cat.a, by = "uniqueID")
SEM.a.df.cat$Disturbance <- as.factor(SEM.b.df.cat$Disturbance)

m2.fit <- sem(m1, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
summary(m2.fit, stand=TRUE, rsq=TRUE)
standardizedSolution(m2.fit, type="std.all")
#Save data
saveRDS(standardizedSolution(m2.fit, type="std.all"), 
        file = here::here("data/SEM_posttransient.rds")) 
  #To read: object <- readRDS(here("data/SEM_posttransient.rds"))

# m2psem <- psem(
#   lm(TStability ~ TVR + TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.a.df),
#   lm(TVR ~ TRichness + TEvenness + Nitrogen + Disturbance + field, data = SEM.a.df),
#   lm(TRichness ~  Nitrogen + Disturbance + field, data = SEM.a.df),
#   lm(TEvenness ~  TRichness + Nitrogen + Disturbance + field, data = SEM.a.df)
# )
# #Path coefficient and significance summary
# summary(m2psem) 
# 
# # bootstrapped model fit for post-transient model
# m2psem_boot <- semEff::bootEff(m2psem, R = 4000, parallel = "multicore", ncpus = 4)
# m2semeff <- semEff::semEff(m2psem_boot)
# #Summary of effects table
# summary(m2semeff)
# 
# #Save bootstrapped fits----
# saveRDS(
#   list(
#     first7 = m1psem_boot,
#     last7 = m2psem_boot
#   ),
#   file = here::here("data/bootstrapped_semfits.rds")
# )

############################
###(In)Direct Effects and Mediation###
############################
#Create a lavaan model with direct effects as before, but add in indirect effects as well

m.indirect <- '#Direct effects on Stab
                TStability ~ c1*Nitrogen + c2*Disturbance+ b2*TEvenness + b4*TRichness + 
                              d1*TVR
               #Mediators
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

#Fit before years model
m.indirect.fit.b <- sem(m.indirect, data=SEM.b.df.cat, se="bootstrap", test="bootstrap")
summary(m.indirect.fit, stand=TRUE, rsq=TRUE) #look at rsq values
standardizedSolution(m.indirect.fit.b, type="std.all")
#Save data
saveRDS(standardizedSolution(m.indirect.fit.b, type="std.all"), 
        file = here::here("data/SEM_indirect_transient.rds")) 
#Fit after years model
m.indirect.fit.a <- sem(m.indirect, data=SEM.a.df.cat, se="bootstrap", test="bootstrap")
summary(m.indirect.fit.a, stand=TRUE, rsq=TRUE)
standardizedSolution(m.indirect.fit.a, type="std.all")
#Save data
saveRDS(standardizedSolution(m.indirect.fit.a, type="std.all"), 
        file = here::here("data/SEM_indirect_posttransient.rds")) 

#Everything below this point was graphing code for old indirect effects plot. With new code,
  #indirect effects will now be individually calculated for select pathways
# #Create theme and plotting functions
# Figtheme<-theme(            axis.line = element_line(colour = "black"),
#                             panel.border = element_rect(colour = NA,fill=NA),
#                             panel.background = element_rect(fill = "white"),
#                             plot.title = element_text(size = 12, hjust = .5),
#                             axis.text = element_text(size = 12, color = "black"),
#                             axis.title = element_text(size = 12),
#                             legend.title = element_blank(),
#                             text = element_text(size = 12))
# number_ticks <- function(n) {function(limits) pretty(limits, n)}
# 
# #############
# 
# #Create dataframe pulling effect sizes from list - Stability before direct
# Stability1 <- data.frame(c(m1semeff[["Effects"]][["TStability"]][["Direct"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TStability"]][["Direct"]][["Disturbance"]], 
#                            m1semeff[["Effects"]][["TStability"]][["Direct"]][["TVR"]],
#                            m1semeff[["Effects"]][["TStability"]][["Direct"]][["TRichness"]],
#                            m1semeff[["Effects"]][["TStability"]][["Direct"]][["TEvenness"]]),
#                          c("Direct","Direct","Direct","Direct","Direct"))%>% 
#   dplyr::rename(value = 1, path = 2) %>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance", "Synchrony", "Richness","Evenness"))
# 
# #Stability before indirect
# Stability2 <- data.frame(c(m1semeff[["Effects"]][["TStability"]][["Indirect"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TStability"]][["Indirect"]][["Disturbance"]],
#                            m1semeff[["Effects"]][["TStability"]][["Indirect"]][["TRichness"]],
#                            m1semeff[["Effects"]][["TStability"]][["Indirect"]][["TEvenness"]]),
#                          c("Indirect","Indirect","Indirect","Indirect")) %>% 
#   dplyr::rename(value = 1, path = 2)%>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance", "Richness","Evenness"))
# 
# Stability <- rbind(Stability1, Stability2) 
# 
# #Stability plot
# p1 <- ggplot(Stability, aes(x=effect,y=value, fill=path)) +
#   geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0), #Set each bar to same width, center around tick
#            width = 0.6) +
#   xlab("") +
#   ylab("Standardized Effect Size") +
#   ggtitle("Effects on Community Stability") +
#   scale_y_continuous(breaks = number_ticks(5)) +
#   scale_fill_manual(values = c("Direct" = "blue", "Indirect" = "purple")) +
#   geom_hline(yintercept=0) +
#   coord_flip() +
#   Figtheme
# 
# ###########
# 
# Synchrony1 <- data.frame(c(m1semeff[["Effects"]][["TVR"]][["Direct"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TVR"]][["Direct"]][["Disturbance"]],
#                            m1semeff[["Effects"]][["TVR"]][["Direct"]][["TRichness"]],
#                            m1semeff[["Effects"]][["TVR"]][["Direct"]][["TEvenness"]]),
#                          c("Direct","Direct","Direct","Direct"))%>% 
#   dplyr::rename(value = 1, path = 2) %>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance", "Richness","Evenness"))
# 
# Synchrony2 <- data.frame(c(m1semeff[["Effects"]][["TVR"]][["Indirect"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TVR"]][["Indirect"]][["Disturbance"]],
#                            m1semeff[["Effects"]][["TVR"]][["Indirect"]][["TRichness"]]),
#                          c("Indirect","Indirect","Indirect")) %>% 
#   dplyr::rename(value = 1, path = 2)%>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance", "Richness"))
# 
# Synchrony <- rbind(Synchrony1, Synchrony2) 
# 
# #Synchrony plot
# p2 <- ggplot(Synchrony, aes(x=effect,y=value, fill=path)) +
#   geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0), #Set each bar to same width, center around tick
#            width = 0.6) +
#   xlab("") +
#   ylab("Standardized Effect Size") +
#   scale_y_continuous(breaks = number_ticks(5)) +
#   ggtitle("Effects on Community Synchrony") +
#   scale_fill_manual(values = c("Direct" = "blue", "Indirect" = "purple")) +
#   geom_hline(yintercept=0) +
#   coord_flip() +
#   Figtheme
# 
# ###########
# 
# Richness1 <- data.frame(c(m1semeff[["Effects"]][["TRichness"]][["Direct"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TVR"]][["Direct"]][["Disturbance"]]),
#                          c("Direct","Direct"))%>% 
#   dplyr::rename(value = 1, path = 2) %>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance"))
# 
# #Richness plot
# p3 <- ggplot(Richness1, aes(x=effect,y=value, fill=path)) +
#   geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0), #Set each bar to same width, center around tick
#            width = 0.15) +
#   xlab("") +
#   ylab("Standardized Effect Size") +
#   scale_y_continuous(breaks = number_ticks(5)) +
#   ggtitle("Effects on Species Richness") +
#   scale_fill_manual(values = c("Direct" = "blue", "Indirect" = "purple")) +
#   geom_hline(yintercept=0) +
#   coord_flip() +
#   Figtheme
# 
# ############
# 
# Evenness1 <- data.frame(c(m1semeff[["Effects"]][["TEvenness"]][["Direct"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TEvenness"]][["Direct"]][["Disturbance"]],
#                            m1semeff[["Effects"]][["TEvenness"]][["Direct"]][["TRichness"]]),
#                          c("Direct","Direct","Direct"))%>% 
#   dplyr::rename(value = 1, path = 2) %>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance","Richness"))
# 
# Evenness2 <- data.frame(c(m1semeff[["Effects"]][["TEvenness"]][["Indirect"]][["Nitrogen"]],
#                            m1semeff[["Effects"]][["TEvenness"]][["Indirect"]][["Disturbance"]]),
#                          c("Indirect","Indirect")) %>% 
#   dplyr::rename(value = 1, path = 2)%>% 
#   dplyr::mutate(effect = c("Nitrogen", "Disturbance"))
# 
# Evenness <- rbind(Evenness1, Evenness2) 
# 
# #Synchrony plot
# p4 <- ggplot(Evenness, aes(x=effect,y=value, fill=path)) +
#   geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0), #Set each bar to same width, center around tick
#            width = 0.6) +
#   xlab("") +
#   ylab("Standardized Effect Size") +
#   scale_y_continuous(breaks = number_ticks(5)) +
#   ggtitle("Effects on Species Evenness") +
#   scale_fill_manual(values = c("Direct" = "blue", "Indirect" = "purple")) +
#   geom_hline(yintercept=0) +
#   coord_flip() +
#   Figtheme
# 
# quartz(width=8, height=5)
# 
# ggarrange(
#   p1, p2, p3, p4, labels = c("A", "B", "C", "D"),
#   common.legend = TRUE, legend = "bottom"
# )

