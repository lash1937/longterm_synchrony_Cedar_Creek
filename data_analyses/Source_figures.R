
################################
# DATA WRANGLING FOR MAIN
# FIGURES 1 - 3 & SUPPLEMENTAL
################################

library(plyr)
library(dplyr)
library(ggplot2)
library(ape)
library(devtools)
library(tidyr)
library(vegan)
library(viridis)
library(stringr)
library(here)
library(codyn)
library(tidyverse)
library(tsvr)
library(ggpubr)
library(investr)
library(ggeffects)
library(patchwork)

source(here::here("data_cleaning/subsetting_CC.R"))

aictable<-function(X,m){
  rnames<-row.names(X)
  AICc<-X$AIC+2*X$df*(X$df+1)/(m-X$df-1)     #small-sample correction
  logL<-X$df-X$AIC/2                         #Log-likelihood
  tab<-data.frame(X[,1],logL,AICc)           #Remove AIC column; add logL and AICc
  colnames(tab)[1]<-c("Params")              #Rename "df" column
  row.names(tab)<-rnames
  tab<-tab[order(tab$AICc),]                 #Sort by ascending AICc value
  deltaAICc<-tab$AICc-min(tab$AICc)          #Delta AICc
  weight<-exp(-deltaAICc/2)/sum(exp(-deltaAICc/2))  #Weights
  cumwt<-weight                              #Column for cumulative weight
  for(i in 2:dim(X)[1]){
    cumwt[i]<-cumwt[i-1]+cumwt[i]              #Accumulate weight from the top
  }
  tab<-data.frame(tab,deltaAICc,weight,cumwt)
  tab<-round(tab,4)
  tab
}
## start with file da.widesynch###

###########################
# Figure 1 data cleaning
###########################

# want data formatted with species, year, location, abundance

unique_ID_exp12 <- tidyr::unite(da.widesynch, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

#converting from wide to longform, keeps the zeros##

unique_ID_long <- tidyr::pivot_longer(unique_ID_exp12, 8:length(unique_ID_exp12),
                                      names_to = "Species", values_to = "Abundance")

unique_ID_long$year<-as.numeric(as.character(unique_ID_long$year))

###just plot information#

####CALCULATING STABILITY#########

st_all <- codyn::community_stability(unique_ID_long, time.var = "year", abundance.var = "Abundance", replicate.var = "uniqueID")

st_all_sep <- tidyr::separate(st_all, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)

st_all_cont<-st_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
st_all_cont$Nitrogen<-mapvalues(st_all_cont$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                                to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))
st_all_cont$Nitrogen<-as.numeric(as.character(st_all_cont$Nitrogen))


##
st_all_cont_minus9 <- st_all_cont%>%     #####dataset for fig 1 B ##### ##got rid of treatment 9 because it doesn't have micronut###
  dplyr::filter(ntrt < 9)

## BIOMASS DATAFRAME ##

biomass_df <- unique_ID_long %>%
  dplyr::group_by(uniqueID, year)%>%
  dplyr::summarize(total_biomass = sum(Abundance))

biomass_overtime <- biomass_df %>%
  dplyr::group_by(uniqueID)%>%
  dplyr::summarise(mean_biomass = mean(total_biomass), stdev_biomass = sd(total_biomass))%>%
  dplyr::mutate(hand_stab = mean_biomass/stdev_biomass)

biomass_all_sep <- tidyr::separate(biomass_overtime, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)

biomass_all_cont<-biomass_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
biomass_all_cont$Nitrogen<-mapvalues(biomass_all_cont$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                                     to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))
biomass_all_cont$Nitrogen<-as.numeric(as.character(biomass_all_cont$Nitrogen))


##
biomass_all_cont_minus9 <- biomass_all_cont%>%     #####dataset for fig 1 B ##### ##got rid of treatment 9 because it doesn't have micronut###
  dplyr::filter(ntrt < 9)



test_stab <- dplyr::left_join(biomass_all_cont_minus9, st_all_cont_minus9)

####CALCULATE VR USING CODYN PACKAGE##########

# let's calculate the variance ratio
# setting average.replicates to false allows us to
VR_all <- codyn::variance_ratio(unique_ID_long, time.var = "year",
                                species.var = "Species", abundance.var = "Abundance",
                                bootnumber = 1, replicate = "uniqueID",
                                average.replicates = FALSE)



VR_all_sep <- tidyr::separate(VR_all, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_")


VR_all_cont<-VR_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
VR_all_cont$Nitrogen<-mapvalues(VR_all_cont$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                                to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))
VR_all_cont$Nitrogen<-as.numeric(as.character(VR_all_cont$Nitrogen))

VR_all_cont_minus9 <- VR_all_cont%>% 
  dplyr::filter(ntrt < 9)     #####dataset for fig 1 A ##### ##got rid of treatment 9 because it doesn't have micronut###
VR_all_cont_minus9 <- VR_all_cont_minus9 %>% 
  tidyr::unite("uniqueID", 1:5, sep="_",remove = FALSE)

VR_all_cont_minus9$disk <- as.factor(VR_all_cont_minus9$disk)
### confidence intervals for graph, disturbed and undisturbed seperate for two linear models ##

#### compare synchrony and stability for the two control treatments ###
VR_all_cont_control <- VR_all_cont %>% 
  dplyr::filter(ntrt %in% c(1,9))

st_all_cont_control <- st_all_cont %>% 
  dplyr::filter(ntrt %in% c(1,9))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#### confidence intervals for the VR ###

VR_all_cont_minus9_E001<-subset(VR_all_cont_minus9, exp==1)

VR_all_cont_minus9_E002<-subset(VR_all_cont_minus9, exp==2)

synchE001lm<-lm(VR~Nitrogen, data = VR_all_cont_minus9_E001)

synchE002lm<-lm(VR~Nitrogen, data = VR_all_cont_minus9_E002)

NitrogenX <- data.frame(
  Nitrogen= seq(0, 30)
)

predE001VR<-as.data.frame(predFit(synchE001lm, newdata = NitrogenX, interval = "confidence", level=0.95))

predE002VR<-as.data.frame(predFit(synchE002lm, newdata = NitrogenX, interval = "confidence", level=0.95))

predE001VR_1<-cbind(NitrogenX, predE001VR)

predE001VR_1$disk<-0

predE002VR_2<-cbind(NitrogenX, predE002VR)

predE002VR_2$disk<-1

confdfVR<-rbind(predE001VR_1, predE002VR_2)

#### confidence intervals for the Stability plot ###

st_all_cont_minus9_E001<-subset(st_all_cont_minus9, exp==1)

st_all_cont_minus9_E002<-subset(st_all_cont_minus9, exp==2)

stabilityE001lm<-lm(stability~Nitrogen, data = st_all_cont_minus9_E001)

stabilityE002lm<-lm(stability~Nitrogen, data = st_all_cont_minus9_E002)

NitrogenX <- data.frame(
  Nitrogen= seq(0, 30)
)

predE001Stability<-as.data.frame(predFit(stabilityE001lm, newdata = NitrogenX, interval = "confidence", level=0.95))

predE002Stability<-as.data.frame(predFit(stabilityE002lm, newdata = NitrogenX, interval = "confidence", level=0.95))

predE001Stability_1<-cbind(NitrogenX, predE001Stability)

predE001Stability_1$disk<-0

predE002Stability_1<-cbind(NitrogenX, predE002Stability)

predE002Stability_1$disk<-1

confdfStability<-rbind(predE001Stability_1, predE002Stability_1)

#######################
# Fig 2
# Time series figure
########################


##unique_ID_long

unique_ID_long2 <- unique_ID_long %>%  tidyr::separate("uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_") %>%  separate("expntrtyear", c("exp2", "ntrt2", "year"), sep="_") %>% select(field, exp, ntrt, plot, year, Species, Abundance )

##remove extra text##

unique_ID_long2$species<-str_remove(unique_ID_long2$Species, "mass.above.")

## summarize biomass all to get most abundance spp ###, 


biomassbysp<-unique_ID_long2 %>% dplyr::group_by(species) %>% dplyr::summarise(totalbiomass=sum(Abundance, na.rm=TRUE))

biomassbysp2 <- merge(biomassbysp,sp.df, by='species', all.x=TRUE)



##pull out the top 5 most abundance species ##


top5sp<-subset(unique_ID_long2, species=="AGROPYRON REPENS"|Species=="POA PRATENSIS"|Species=="ARTEMISIA LUDOVICIANA"|Species=="SCHIZACHYRIUM SCOPARIUM"|Species=="SOLIDAGO RIGIDA")

topbyfunc<-subset(unique_ID_long2, species=="AGROPYRON REPENS"|species=="POA PRATENSIS"|species=="ARTEMISIA LUDOVICIANA"|species=="SCHIZACHYRIUM SCOPARIUM"|species=="SOLIDAGO RIGIDA"|species=="RUBUS SP."|species=="LATHYRUS VENOSUS"|species=="POLYGONUM CONVOLVULUS")


totalbiomass<-unique_ID_long2 %>% dplyr::group_by(field, exp,year, plot, ntrt) %>% dplyr::summarise(totalbiomass=sum(Abundance, na.rm=TRUE))


##average across treatments##

top5sp_avg<-top5sp%>% dplyr::group_by(exp,year,ntrt, species) %>% dplyr::summarise(meanbiomass=mean(Abundance))

topbyfunc_avg<-topbyfunc%>% dplyr::group_by(exp,year,ntrt, species) %>% dplyr::summarise(meanbiomass=mean(Abundance))


totalbiomass_avg<-totalbiomass %>% dplyr::group_by(exp,year,ntrt) %>% dplyr::summarise(meantotalbiomass=mean(totalbiomass))


##subset to only control and high nutrient###

top5sp_avg_subset<-subset(top5sp_avg, ntrt==1|ntrt==6)

topbyfunc_subset<-subset<-subset(topbyfunc_avg, ntrt==1|ntrt==6)

top5sp_avg_subset$ntrt <- factor(top5sp_avg_subset$ntrt, levels = c("1", "6"))

topbyfunc_subset$ntrt <- factor(topbyfunc_subset$ntrt, levels = c("1", "6"))


totalbiomass_subset<-subset(totalbiomass_avg, ntrt==1|ntrt==6)

totalbiomass_subset$ntrt <- factor(totalbiomass_subset$ntrt, levels = c("1", "6"))

totalbiomass_subset<-totalbiomass_subset %>% dplyr::mutate(species="total biomass")






##############
# FIGURE 3
##############

## sub-setting synchrony calculations into exp 1 and exp 2 to handle missing years in exp 2

###E001 #####

sub_exp1_long <- da.widesynch %>% 
  dplyr::filter(exp == "1") %>% 
  tidyr::unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = FALSE)%>% 
  dplyr::select(-subplot, -ntrt2, -expntrtyear, -expntrtfieldyear)

sub_exp1_long <- sub_exp1_long %>% 
  tidyr::pivot_longer(9:length(sub_exp1_long),
                      names_to = "Species", values_to = "Abundance")%>%filter(Abundance > 0)

sub_exp1_wide <- tidyr::pivot_wider(sub_exp1_long, names_from = year, values_from = Abundance)%>%select(-Species,-field,-exp,-plot,-disk,-ntrt,-expfieldplot)

sub_exp1_wide[is.na(sub_exp1_wide)] <- 0

plotnames<-c(unique(sub_exp1_long$uniqueID))

plotnames2<-str_sort(plotnames)

vreq_exp_1<-sub_exp1_wide

matrix3<-matrix(NA, ncol=4, nrow=length(plotnames))   ### creating a blank matrix####
colnames(matrix3)<-c("plot", "classicVR", "comm", "pop")
plotnames_exp1<-c(unique(vreq_exp_1$uniqueID))

##### Loop is going through each unique ID (plot) and calculating pop variability, comm variability, and total variance ratio#####

for(s in 1:length(plotnames_exp1)){
  current.plot<-plotnames_exp1[s]
  temp<-vreq_exp_1 %>% 
    dplyr::filter(uniqueID %in% current.plot) %>% 
    dplyr::select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  matrix3[s,2] <- vr[[3]]
  matrix3[s,1] <- current.plot
  matrix3[s,3]<- vr[[1]]
  matrix3[s,4]<- vr[[2]]
}


popcomvar_exp_1 <- as.data.frame(matrix3)
####dataframe for exp 1###

### E002####
sub_exp2_long <- da.widesynch %>% 
  dplyr::filter(exp == "2") %>% 
  tidyr::unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = FALSE)%>% 
  dplyr::select(-subplot, -ntrt2, -expntrtyear, -expntrtfieldyear)

sub_exp2_long <- sub_exp2_long %>% tidyr::pivot_longer(9:length(sub_exp2_long),
                                                       names_to = "Species", values_to = "Abundance")%>%filter(Abundance > 0)

sub_exp2_wide <- tidyr::pivot_wider(sub_exp2_long, names_from = year, values_from = Abundance)%>%select(-Species,-field,-exp,-plot,-disk,-ntrt,-expfieldplot)

sub_exp2_wide[is.na(sub_exp2_wide)] <- 0

plotnames<-c(unique(sub_exp2_long$uniqueID))

plotnames2<-str_sort(plotnames)

matrix1<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix1)<-c("plot", "classic", "short", "long")

vreq_exp_2<-sub_exp2_wide

matrix4<-matrix(NA, ncol=4, nrow=length(plotnames))   ### creating a blank matrix####
colnames(matrix4)<-c("plot", "classicVR", "comm", "pop") # add mean and stdev columns
plotnames_exp2<-c(unique(vreq_exp_2$uniqueID))


for(s in 1:length(plotnames_exp2)){
  current.plot<-plotnames_exp2[s]
  temp<-vreq_exp_2 %>% 
    dplyr::filter(uniqueID %in% current.plot) %>% 
    dplyr::select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  matrix4[s,2] <- vr[[3]]
  matrix4[s,1] <- current.plot
  matrix4[s,3]<- vr[[1]]
  matrix4[s,4]<- vr[[2]]
}


popcomvar_exp_2 <- as.data.frame(matrix4)

# combine E001 and E002 dataframes###

popcomvar_exp12<-rbind(popcomvar_exp_1, popcomvar_exp_2)

popcomvar_exp12_2<- tidyr::separate(popcomvar_exp12, "plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

popcomvar_exp12_2<-popcomvar_exp12_2 %>% dplyr::mutate(Nitrogen=ntrt)
popcomvar_exp12_2$Nitrogen<-mapvalues(popcomvar_exp12_2$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                                      to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))
popcomvar_exp12_2$Nitrogen<-as.numeric(as.character(popcomvar_exp12_2$Nitrogen))


### average across the disk and nutrients####
popcomvar_exp12_2$disk<-as.factor(popcomvar_exp12_2$disk)
popcomvar_exp12_2$Nitrogen<-as.factor(popcomvar_exp12_2$Nitrogen)
popcomvar_exp12_2$classicVR<-as.numeric(popcomvar_exp12_2$classicVR)
popcomvar_exp12_2$comm<-as.numeric(popcomvar_exp12_2$comm)
popcomvar_exp12_2$pop<-as.numeric(popcomvar_exp12_2$pop)

# summarize (mean_of_stab = mean(mean), std_of_stab=mean(std))
avgpopcommvar<- popcomvar_exp12_2 %>% dplyr::group_by(disk,Nitrogen) %>% 
  dplyr::summarize(meanVR = mean(classicVR), meanpop=mean(pop), meancomm=mean(comm))

# remove outliers for plotting purposes
biomass_all_cont_minus9andoutliers <- biomass_all_cont_minus9

avg_test_stab <- test_stab %>%
  dplyr::group_by(disk, Nitrogen)%>%
  dplyr::summarise(mean_st = mean(stability))

# still calculate averages with outliers included
biomass_disk_N <- biomass_all_cont_minus9 %>%
  dplyr::group_by(disk, Nitrogen)%>%
  dplyr::summarise(meanofmean_biomass = mean(mean_biomass), meanofsd_biomass = mean(stdev_biomass))%>%
  dplyr::mutate(hand_stab = meanofmean_biomass/meanofsd_biomass)

avg_test <- dplyr::left_join(avg_test_stab, biomass_disk_N)


biomass_disk_N$disk<-as.factor(biomass_disk_N$disk)
biomass_disk_N$Nitrogen<-as.factor(biomass_disk_N$Nitrogen)
biomass_all_cont_minus9andoutliers$disk<-as.factor(biomass_all_cont_minus9andoutliers$disk)
biomass_all_cont_minus9andoutliers$Nitrogen<-as.factor(biomass_all_cont_minus9andoutliers$Nitrogen)




#