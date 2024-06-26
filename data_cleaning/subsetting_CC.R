#####################################
# This script subsets and tidies data from the Cedar Creek E001 and E002 
# datasets, for use in analyses and figure creation. 
# Four resulting datasets are created: 1. before_years, which is a wide 
# dataframe of all species-level info from the first 7 years, 2. SEM.b.df, 
# which compiles calculations of synchrony, stability, richness, and evenness 
# for each unique plot from the first 7 years, 3. after_years, which is a wide 
# dataframe of all species-level info from the last 7 years with data(1994, 1996, 
# 1997, 1999, 2000, 2002, 2004), and 4. SEM.a.df, which compiles calculations 
# of synchrony, stability, richness, and evenness for each unique plot from 
# the last 7 years.
#####################################

library(here)
library(codyn)

#Read in data and functions
source(here::here("data_cleaning/cleaning_CC.R"))
source(here::here("Functions/functions.R"))

###All Years Data###
#Sum total biomass of each species per plot over all the years
unique_ID_exp12 <- tidyr::unite(da.widesynch, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

unique_ID_long <- tidyr::pivot_longer(unique_ID_exp12, 8:length(unique_ID_exp12),
                               names_to = "Species", values_to = "Abundance")

unique_ID_long$Abundance<-as.numeric(unique_ID_long$Abundance)
unique_ID_long$year<-as.numeric(unique_ID_long$year)

#Diversity-----
Richness.df <- unique_ID_long %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Rich = sum(Abundance > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Richness = mean(Rich))

Evenness.df <- unique_ID_long %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Even = Evar(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Evenness = mean(Even, na.rm=T))

#Synchrony-----
VR_all <- variance_ratio(unique_ID_long, time.var = "year",
                         species.var = "Species", abundance.var = "Abundance",
                         bootnumber = 1, replicate = "uniqueID",
                         average.replicates = FALSE)

#Stability-----
stab_all <- community_stability(unique_ID_long, time.var = "year",
                                abundance.var = "Abundance", replicate.var = "uniqueID")

#All Variables------
#Create consolidated dataframe of diversity metrics, VR, stability, disturbance, and nutrients
SEM.df <- dplyr::left_join(Richness.df, Evenness.df)%>%
  dplyr::left_join(VR_all)%>%
  dplyr::left_join(stab_all)%>%
  dplyr::select(uniqueID, Richness, Evenness, VR, stability)%>%
  tidyr::separate(uniqueID, into = c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  dplyr::select(-c(plot)) %>% 
  dplyr::rename(Disturbance = disk, Nutrients = ntrt, Stability = stability) %>% 
  mutate(grid = factor(paste0(field, exp))) %>% 
  dplyr::relocate(grid, .after = exp)

#Add in dummy variables for Field effects
field <- as.data.frame(model.matrix(~ field, data = SEM.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.df$uniqueID)

SEM.df <- left_join(SEM.df, field, by = "uniqueID")


#Convert Nutrients to Continuous Values-----
SEM.df <- SEM.df %>% dplyr::mutate(Micronut=Nutrients, Nitrogen=Nutrients)
SEM.df$Micronut<-mapvalues(SEM.df$Micronut, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"), 
                           to=c("1", "1", "1" ,"1", "1", "1", "1", "1", "0"))
SEM.df$Nitrogen<-mapvalues(SEM.df$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                           to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))
SEM.df$Nitrogen<-as.numeric(as.character(SEM.df$Nitrogen))
SEM.df$Disturbance <- as.factor(SEM.df$Disturbance)

#Remove Nutrient Control (Micronut = 0)
SEM.df <- subset(SEM.df, Micronut!=0)


#Transient Years Data Subset-----
before_years <- subset(da.widesynch, as.numeric(year) <= 7)
before_years$years <- droplevels(before_years$year)

##Transient years dataframe-----
unique_ID_exp12_before <- tidyr::unite(before_years, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

unique_ID_long_before <- tidyr::pivot_longer(unique_ID_exp12_before, cols = 8:(length(unique_ID_exp12_before)-1),
                                      names_to = "Species", values_to = "Abundance")

unique_ID_long_before$Abundance<-as.numeric(unique_ID_long_before$Abundance)
unique_ID_long_before$year<-as.numeric(unique_ID_long_before$year)

#Diversity-----
Richness.b.df <- unique_ID_long_before %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Rich = sum(Abundance > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Richness = mean(Rich))

Evenness.b.df <- unique_ID_long_before %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Even = Evar(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Evenness = mean(Even, na.rm=T))

#Synchrony-----
VR_before <- variance_ratio(unique_ID_long_before, time.var = "year",
                            species.var = "Species", abundance.var = "Abundance",
                            bootnumber = 1, replicate = "uniqueID",
                            average.replicates = FALSE)

#Population & Community Variability-----
sub_bf_wide <- unique_ID_long_before %>% 
  dplyr::select(uniqueID, years, Species, Abundance) %>% 
  tidyr::pivot_wider(names_from = years, values_from = Abundance) %>%
  dplyr::select(-Species)

plotnames <- c(unique(sub_bf_wide$uniqueID))
matrix_bf <- matrix(NA, ncol=4, nrow=length(plotnames))  
colnames(matrix_bf) <- c("uniqueID", "classicVR", "comm", "pop")

for(s in 1:length(plotnames)){
  current.plot <- plotnames[s]
  temp <- sub_bf_wide %>% 
    dplyr::filter(uniqueID %in% current.plot) %>% 
    dplyr::select(-uniqueID)
  temp<-as.matrix(temp)
  temp <- temp[rowSums(temp[])>0,]
  vr<-vreq_classic(temp)
  matrix_bf[s,2] <- vr[[3]]
  matrix_bf[s,1] <- current.plot
  matrix_bf[s,3]<- vr[[1]]
  matrix_bf[s,4]<- vr[[2]]
}
popcomvar_bf <- as.data.frame(matrix_bf)

#Stability-----
stab_before <- community_stability(unique_ID_long_before, time.var = "year",
                                   abundance.var = "Abundance", replicate.var = "uniqueID")

#Total Biomass-----
biomass_overtime_bf <- unique_ID_long_before %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(total_biomass = sum(Abundance)) %>% 
  dplyr::group_by(uniqueID)%>%
  dplyr::summarise(mean_biomass = mean(total_biomass), stdev_biomass = sd(total_biomass))%>%
  dplyr::mutate(hand_stab = mean_biomass/stdev_biomass)

#All Transient Variables-----
# Create consolidated dataframe of diversity metrics, VR, stability, population and 
# community variance, mean biomass, disturbance, and nutrients

SEM.b.df <- left_join(Richness.b.df, Evenness.b.df) %>%
  dplyr::left_join(VR_before)%>%
  dplyr::left_join(stab_before)%>%
  dplyr::left_join(biomass_overtime_bf)%>%
  dplyr::left_join(popcomvar_bf) %>% 
  dplyr::select(uniqueID, Richness, Evenness, VR, stability, mean_biomass, stdev_biomass,
                comm, pop)%>%
  tidyr::separate(uniqueID, into = c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  dplyr::select(-c(plot)) %>% 
  dplyr::rename(Disturbance = disk, Nutrients = ntrt, Stability = stability) %>% 
  mutate(grid = factor(paste0(field, exp))) %>% 
  dplyr::relocate(grid, .after = exp)
  
#Add in dummy variables for Field effects
field <- as.data.frame(model.matrix(~ field, data = SEM.b.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.b.df$uniqueID)

SEM.b.df <- left_join(SEM.b.df, field, by = "uniqueID")
  

#Convert Nutrients to Continuous Values-----
SEM.b.df <- SEM.b.df %>% dplyr::mutate(Micronut=Nutrients, Nitrogen=Nutrients)
SEM.b.df$Micronut<-mapvalues(SEM.b.df$Micronut, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"), 
                             to=c("1", "1", "1" ,"1", "1", "1", "1", "1", "0"))
SEM.b.df$Nitrogen<-mapvalues(SEM.b.df$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                             to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))

SEM.b.df$Nitrogen <- as.numeric(as.character(SEM.b.df$Nitrogen))
SEM.b.df$pop <- as.numeric(SEM.b.df$pop)
SEM.b.df$comm <- as.numeric(SEM.b.df$comm)
SEM.b.df$Disturbance <- as.numeric(SEM.b.df$Disturbance)

#Remove Nutrient Control (Micronut = 0)
SEM.b.df <- SEM.b.df %>%
  subset(Micronut!=0)


#Post-transient Years Data Subset-----
after_years <- subset(da.widesynch, year %in% c(1994, 1996, 1997, 1999, 2000, 2002, 2004))
after_years$years <- droplevels(after_years$year)

#Post-transient years dataframe-----
unique_ID_exp12_after <- tidyr::unite(after_years, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

unique_ID_long_after <- tidyr::pivot_longer(unique_ID_exp12_after, 8:(length(unique_ID_exp12_after)-1),
                                     names_to = "Species", values_to = "Abundance")

unique_ID_long_after$Abundance<-as.numeric(unique_ID_long_after$Abundance)
unique_ID_long_after$year<-as.numeric(unique_ID_long_after$year)

#Diversity-----
Richness.a.df <- unique_ID_long_after %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Rich = sum(Abundance > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Richness = mean(Rich))

Evenness.a.df <- unique_ID_long_after %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Even = Evar(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Evenness = mean(Even, na.rm=T))

#Synchrony-----
VR_after <- variance_ratio(unique_ID_long_after, time.var = "year",
                           species.var = "Species", abundance.var = "Abundance",
                           bootnumber = 1, replicate = "uniqueID",
                           average.replicates = FALSE)

#Population & Community Variability-----
sub_af_wide <- unique_ID_long_after %>% 
  dplyr::select(uniqueID, years, Species, Abundance) %>% 
  tidyr::pivot_wider(names_from = years, values_from = Abundance)%>%
  dplyr::select(-Species)

plotnames<-c(unique(sub_af_wide$uniqueID))
matrix_af<-matrix(NA, ncol=4, nrow=length(plotnames))  
colnames(matrix_af)<-c("uniqueID", "classicVR", "comm", "pop")

for(s in 1:length(plotnames)){
  current.plot <- plotnames[s]
  temp <- sub_af_wide %>% 
   dplyr:: filter(uniqueID %in% current.plot) %>% 
   dplyr::select(-uniqueID)
  temp <- as.matrix(temp)
  temp <- temp[rowSums(temp[])>0,]
  vr <- vreq_classic(temp)
  matrix_af[s,2] <- vr[[3]]
  matrix_af[s,1] <- current.plot
  matrix_af[s,3]<- vr[[1]]
  matrix_af[s,4]<- vr[[2]]
}

popcomvar_af <- as.data.frame(matrix_af)

#Stability
stab_after <- community_stability(unique_ID_long_after, time.var = "year",
                                  abundance.var = "Abundance", replicate.var = "uniqueID")

#Total Biomass
biomass_overtime_af <- unique_ID_long_after %>%
  dplyr::group_by(uniqueID, year)%>%
  dplyr::summarize(total_biomass = sum(Abundance)) %>% 
  group_by(uniqueID) %>%
  dplyr::summarise(mean_biomass = mean(total_biomass), stdev_biomass = sd(total_biomass)) %>%
  dplyr::mutate(hand_stab = mean_biomass/stdev_biomass)

#All Post-transient Variables
# Create consolidated dataframe of diversity metrics, VR, stability, population and 
# community variance, mean biomass, disturbance, and nutrients

SEM.a.df <- left_join(Richness.a.df, Evenness.a.df) %>%
  dplyr::left_join(VR_after)%>%
  dplyr::left_join(stab_after)%>%
  dplyr::left_join(biomass_overtime_af)%>%
  dplyr::left_join(popcomvar_af) %>% 
  dplyr::select(uniqueID, Richness, Evenness, VR, stability, mean_biomass, stdev_biomass,
                comm, pop)%>%
  tidyr::separate(uniqueID, into = c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  dplyr::select(-c(plot)) %>% 
  dplyr::rename(Disturbance = disk, Nutrients = ntrt, Stability = stability) %>% 
  mutate(grid = factor(paste0(field, exp))) %>% 
  dplyr::relocate(grid, .after = exp)

#Add in dummy variables for Field effects
field.a <- as.data.frame(model.matrix(~ field, data = SEM.a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.a.df$uniqueID)

SEM.a.df <- left_join(SEM.a.df, field.a, by = "uniqueID")

#Convert Nutrients to Continuous Values-----
SEM.a.df<-SEM.a.df %>% dplyr::mutate(Micronut=Nutrients, Nitrogen=Nutrients)
SEM.a.df$Micronut<-mapvalues(SEM.a.df$Micronut, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"), 
                             to=c("1", "1", "1" ,"1", "1", "1", "1", "1", "0"))
SEM.a.df$Nitrogen<-mapvalues(SEM.a.df$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                             to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))

SEM.a.df$Nitrogen<-as.numeric(as.character(SEM.a.df$Nitrogen))
SEM.a.df$pop <- as.numeric(SEM.a.df$pop)
SEM.a.df$comm <- as.numeric(SEM.a.df$comm)
SEM.a.df$Disturbance <- as.numeric(SEM.a.df$Disturbance)

#Remove Nutrient Control (Micronut = 0)
SEM.a.df <- SEM.a.df %>%
  subset(Micronut!=0) 

##################
#Create a dataset with a 10 year window instead of 7 years
##################

#Transient 10 Years Data Subset-----
before_years <- subset(da.widesynch, as.numeric(year) <= 10)
before_years$years <- droplevels(before_years$year)

##Transient 10 years dataframe-----
unique_ID_exp12_10before <- tidyr::unite(before_years, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

unique_ID_long_10before <- tidyr::pivot_longer(unique_ID_exp12_10before, cols = 8:(length(unique_ID_exp12_10before)-1),
                                             names_to = "Species", values_to = "Abundance")

unique_ID_long_10before$Abundance<-as.numeric(unique_ID_long_10before$Abundance)
unique_ID_long_10before$year<-as.numeric(unique_ID_long_10before$year)

#Diversity-----
Richness.10b.df <- unique_ID_long_10before %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Rich = sum(Abundance > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Richness = mean(Rich))

Evenness.10b.df <- unique_ID_long_10before %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Even = Evar(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Evenness = mean(Even, na.rm=T))

#Synchrony-----
VR_10before <- variance_ratio(unique_ID_long_10before, time.var = "year",
                            species.var = "Species", abundance.var = "Abundance",
                            bootnumber = 1, replicate = "uniqueID",
                            average.replicates = FALSE)

#Stability-----
stab_10before <- community_stability(unique_ID_long_10before, time.var = "year",
                                   abundance.var = "Abundance", replicate.var = "uniqueID")


#All Transient 10 year window Variables-----
# Create consolidated dataframe of diversity metrics, VR, stability, population and 
# community variance, mean biomass, disturbance, and nutrients

SEM.10b.df <- left_join(Richness.10b.df, Evenness.10b.df) %>%
  dplyr::left_join(VR_10before)%>%
  dplyr::left_join(stab_10before)%>%
  dplyr::select(uniqueID, Richness, Evenness, VR, stability)%>%
  tidyr::separate(uniqueID, into = c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  dplyr::select(-c(plot)) %>% 
  dplyr::rename(Disturbance = disk, Nutrients = ntrt, Stability = stability) %>% 
  mutate(grid = factor(paste0(field, exp))) %>% 
  dplyr::relocate(grid, .after = exp)

#Add in dummy variables for Field effects
field10 <- as.data.frame(model.matrix(~ field, data = SEM.10b.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.10b.df$uniqueID)

SEM.10b.df <- left_join(SEM.10b.df, field10, by = "uniqueID")


#Convert Nutrients to Continuous Values-----
SEM.10b.df <- SEM.10b.df %>% dplyr::mutate(Micronut=Nutrients, Nitrogen=Nutrients)
SEM.10b.df$Micronut<-mapvalues(SEM.10b.df$Micronut, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"), 
                             to=c("1", "1", "1" ,"1", "1", "1", "1", "1", "0"))
SEM.10b.df$Nitrogen<-mapvalues(SEM.10b.df$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                             to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))

SEM.10b.df$Nitrogen <- as.numeric(as.character(SEM.10b.df$Nitrogen))
SEM.10b.df$Disturbance <- as.numeric(SEM.10b.df$Disturbance)

#Remove Nutrient Control (Micronut = 0)
SEM.10b.df <- SEM.10b.df %>%
  subset(Micronut!=0)



#Post-transient 10 Years Data Subset-----
after_years <- subset(da.widesynch, year %in% c(1991 ,1992, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004))
after_years$years <- droplevels(after_years$year)

#Post-transient years dataframe-----
unique_ID_exp12_10after <- tidyr::unite(after_years, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

unique_ID_long_10after <- tidyr::pivot_longer(unique_ID_exp12_10after, 8:(length(unique_ID_exp12_10after)-1),
                                            names_to = "Species", values_to = "Abundance")

unique_ID_long_10after$Abundance<-as.numeric(unique_ID_long_10after$Abundance)
unique_ID_long_10after$year<-as.numeric(unique_ID_long_10after$year)

#Diversity-----
Richness.10a.df <- unique_ID_long_10after %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Rich = sum(Abundance > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Richness = mean(Rich))

Evenness.10a.df <- unique_ID_long_10after %>%
  dplyr::group_by(uniqueID, year) %>%
  dplyr::summarise(Even = Evar(Abundance)) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(uniqueID) %>%
  dplyr::summarise(Evenness = mean(Even, na.rm=T))

#Synchrony-----
VR_10after <- variance_ratio(unique_ID_long_10after, time.var = "year",
                           species.var = "Species", abundance.var = "Abundance",
                           bootnumber = 1, replicate = "uniqueID",
                           average.replicates = FALSE)

#Stability
stab_10after <- community_stability(unique_ID_long_10after, time.var = "year",
                                  abundance.var = "Abundance", replicate.var = "uniqueID")


#All Post-transient 10 year Variables
# Create consolidated dataframe of diversity metrics, VR, stability, population and 
# community variance, mean biomass, disturbance, and nutrients

SEM.10a.df <- left_join(Richness.10a.df, Evenness.10a.df) %>%
  dplyr::left_join(VR_10after)%>%
  dplyr::left_join(stab_10after)%>%
  dplyr::select(uniqueID, Richness, Evenness, VR, stability)%>%
  tidyr::separate(uniqueID, into = c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  dplyr::select(-c(plot)) %>% 
  dplyr::rename(Disturbance = disk, Nutrients = ntrt, Stability = stability) %>% 
  mutate(grid = factor(paste0(field, exp))) %>% 
  dplyr::relocate(grid, .after = exp)

#Add in dummy variables for Field effects
field.10a <- as.data.frame(model.matrix(~ field, data = SEM.10a.df)) %>% 
  dplyr::rename(fieldA = `(Intercept)`) %>% 
  dplyr::mutate(uniqueID = SEM.10a.df$uniqueID)

SEM.10a.df <- left_join(SEM.10a.df, field.10a, by = "uniqueID")

#Convert Nutrients to Continuous Values-----
SEM.10a.df<-SEM.10a.df %>% dplyr::mutate(Micronut=Nutrients, Nitrogen=Nutrients)
SEM.10a.df$Micronut<-mapvalues(SEM.10a.df$Micronut, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"), 
                             to=c("1", "1", "1" ,"1", "1", "1", "1", "1", "0"))
SEM.10a.df$Nitrogen<-mapvalues(SEM.10a.df$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "7", "8", "9"),
                             to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "17", "27.2", "0.0"))

SEM.10a.df$Nitrogen<-as.numeric(as.character(SEM.10a.df$Nitrogen))
SEM.10a.df$Disturbance <- as.numeric(SEM.10a.df$Disturbance)

#Remove Nutrient Control (Micronut = 0)
SEM.10a.df <- SEM.10a.df %>%
  subset(Micronut!=0) 

