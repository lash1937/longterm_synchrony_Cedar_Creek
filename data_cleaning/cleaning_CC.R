# Reading in and cleaning Cedar Creek data
# both above and belowground for e001 and e002

library(plyr)
library(dplyr)
library(ggplot2)
library(ape)
library(devtools)
library(tidyr)
library(vegan)
library(viridis)
library(stringr)
library(codyn)
library(tidyverse)
library(tsvr)
library(ggpubr)
library(lubridate)
library(here)


# Reading in and cleaning Cedar Creek data###

#Read in e001 data## Plot data for the unplowed fields##

d1s <- read.csv(here::here("data/e001-soil-cn-2019-09-13.csv"), header=T, strip.white=T, skip=1, , fileEncoding = "UTF-8-BOM")

# Delete some extra variables and add dummies to stack exp1 & exp2 data
d1s$natm.nadd <- NULL
d1s$fertadd <- NULL
d1s$burn.trt.1992 <- NA
d1s$ntrt.last.year <- 'CurrentYear'
d1s$ntrt.origin <- 1 # No cessation in E001
d1s$disk <- 0 # E001 was not disked before experiment (disturbed vegetation)
d1s$exp <- 1

#Read in e002 data ##Plot data for the plowed fields##

d2s <- read.csv(here::here("data/e002-soil-cn-2019-09-13.csv"), header=T, strip.white=T, skip=1, fileEncoding = "UTF-8-BOM")

## Delete some extra variables and add dummies to stack exp1 & exp2 data
d2s$X <- NULL
d2s$date <- NULL
d2s$BurnRecBeginning1992 <- NULL
d2s$NtrtRecBefore1992 <- NULL
d2s$TrtRecAfter1991 <- NULL
d2s$TrtRecAfter2013 <- NULL
d2s$burn <- NA
d2s$fence <- NA
d2s$fence.origin <- NA
d2s$depth <- NA
d2s$disk <- 1 # E002 was disked before experiment (disturbed vegetation)
d2s$exp <- 2


#Additioanl plot info for both experiments##

d3s <- read.csv(here::here("data/E001-E002-Soil-CN-2018.csv"), header=T, strip.white=T, skip=0, na.strings=c(""), fileEncoding = "UTF-8-BOM")


# Delete samples that have errors
d3s <- d3s[is.na(d3s$FLAGGED.FOR.RERUN..S.Spillage..E.Error.),]

# Delete some extra variables and add dummies to stack exp1 & exp2 data
d3s$FLAGGED.FOR.RERUN..S.Spillage..E.Error. <- NULL
d3s$Sample <- NULL
d3s$X <- NULL
d3s$burn <- NA
d3s$burn.origin <- NA
d3s$burn.trt.1992 <- NA
d3s$fence <- NA
d3s$fence.origin <- NA
d3s$depth <- NA
d3s$disk <- NA
d3s$nadd <- NA
d3s$ntrt <- NA
d3s$ntrt.last.year <- NA
d3s$ntrt.origin <- NA
d3s$trt.origin <- NA  

d3s$disk[d3s$exp==1] <- 0 # E001 was not disked before experiment (intact vegetation)
d3s$disk[d3s$exp==2] <- 1 # E002 was disked before experiment (disturbed vegetation)

# Check that data sets have the same variables
sort(names(d1s))
sort(names(d2s))
sort(names(d3s))

# Stack exp1 & 2 data
ds <- rbind(d1s,d2s,d3s)


#Read in e001 aboveground plant biomass data ##

d1a <- read.csv(here::here("data/e001-aboveground-mass-2019-09-13.csv"), header=T, strip.white=T, skip=1, fileEncoding = "UTF-8-BOM")



# Delete some extra variables and add dummies to stack exp1 & exp2 data
d1a$natm.nadd <- NULL
d1a$fertadd <- NULL
d1a$burn.trt.1992 <- NA
d1a$ntrt.last.year <- "CurrentYear" # No cessation in E001
d1a$ntrt.origin <- 1 # No cessation in E001
d1a$subplot <- NA
d1a$disk <- 0 # E001 was not disked before experiment (disturbed vegetation)
d1a$exp <- 1


#Read in e002 aboveground plant biomass data ##
d2a <- read.csv(here::here("data/e002-aboveground-mass-2019-09-13.csv"), header=T, strip.white=T, skip=1, fileEncoding = "UTF-8-BOM")


# Delete some extra variables and add dummies to stack exp1 & exp2 data
d2a$X <- NULL
d2a$date <- NULL
d2a$burn <- NA
d2a$fence <- NA
d2a$fence.origin <- NA
d2a$disk <- 1 # E002 was disked before experiment (intact vegetation)
d2a$exp <- 2


# Check that data sets have the same variables
sort(names(d1a))
sort(names(d2a))


# Stack exp1 & 2 data

da <- rbind(d1a,d2a)

# Delete records with key missing data
da <- da[!is.na(da$mass),]

# Disking not done if field D
da <- da[da$field != 'D',]
da$field <- as.factor(as.character(da$field))


# If subplot is missing specify it as a whole plot
da$subplot[is.na(da$subplot)] <- "Whole"

# Code other Nutrient additions
da$other.add <- 1
da$other.add[da$ntrt == 9] <- 0


### da = long form version of full dataset including plots without original treatments###

################################################################################
# Here we make a clean design file of original treatments                      #
# Make clean subplot and year level file with original treatments only.        #
################################################################################


names(da)

# SUBPLOT SCALE Get a list of plots that have original treatments
design.df <- ddply(da, .(field, exp, plot, subplot, ntrt, nadd, disk, other.add, ntrt.origin), colwise(mean, .(mass.above)))
design.df$mass.above <- NULL
dim(design.df)
with(design.df, table(plot, field, exp))
dim(design.df)  ### 603 X 9 ##




# PLOT SCALE Get a list of plots that have original treatments
design.df.plot <- ddply(da, .(field, exp, plot, ntrt, nadd, disk, other.add, ntrt.origin), colwise(mean, .(mass.above)))
design.df.plot$mass.above <- NULL
dim(design.df.plot)
with(design.df.plot, table(plot, field, exp))
dim(design.df.plot)


# YEAR AND SUBPLOT SCALE Get a list of plots that have original treatments
design.df.yr <- ddply(da, .(field, exp, plot, subplot, year, ntrt, disk, other.add, ntrt.origin, burn.origin, fence.origin), colwise(mean, .(mass.above)))
design.df$mass.above <- NULL
dim(design.df)
with(design.df, table(plot, field, exp))
dim(design.df)

#### Important note: burn.origin column isn't exactly right for B E002. The burn doesn't happen until 1992, so should be 54 plots from 1982-1991 and 27 plots 1992-2004###

summary(design.df.yr)
design.orig <- design.df.yr   #dim 9086 X 12

# Delete cessations plots
design.orig <- design.orig[design.orig$ntrt.origin==1,]  #dim #7799   12#


# Delete subplots with experimental burns (after 1992)
#design.orig <- design.orig[design.orig$burn.origin==1,]
design.orig<-subset(design.orig, year < 1992| year >= 1992 & burn.origin==1) #dim #6153   12#


# Fences removed in 2004 (partial removal in field C but still open generally)
# After 2004 all of E002 is unfenced
design.orig$fence[design.orig$year <= 2004 & design.orig$exp==2] <- 1
design.orig$fence[design.orig$year > 2004 & design.orig$exp==2] <- 0

# After 2004 all of E001 is unfenced except in field C
design.orig$fence[design.orig$year <= 2004 & design.orig$exp==1 & design.orig$field != 'C'] <- 1
design.orig$fence[design.orig$year > 2004 & design.orig$exp==2 & design.orig$field != 'C'] <- 0

# Pull out data that is fenced after 2004 as these were invidual plots 
# fenced in field C (I think) after whole field fences were removed. 
design.orig$sel <- TRUE
design.orig$sel[design.orig$year > 2004 & design.orig$fence == 1] <- FALSE
design.orig <- design.orig[design.orig$sel,]
with(design.orig, table(year,fence, exp))
with(design.orig, table(year, field,fence.origin, exp))
with(design.orig, table(year, fence, field, exp))

dim(design.orig)   #dim #6153   12#
summary(design.orig)

# Get field scale burn record
df.burn <- ddply(d1a, .(field, year), colwise(max, .(burn)), na.rm=T)
# Replace original burn record with field summary
da$burn <- NULL
da <- merge(da, df.burn, by=c("field", "year"), all.x=TRUE)


# Merge design.orig file, to extract files for which treatments have not changed
dim(da)  ##72737  X  19### ##includes N cessation and burned plots###
dim(design.orig)  # Unique plot_years with correct treatments### ##dim 6396##

design.orig2<-design.orig %>% dplyr::select(field, year, exp, disk, plot, subplot, ntrt, other.add, ntrt.origin, burn.origin, fence.origin)

da_min<-da %>% dplyr::select(field, year, exp, disk, plot, subplot, ntrt,  species, mass.above) #72737  X 9### ##includes N cessation and burned plots###

da_orig<-merge(design.orig2, da_min, by =c("field", "year", "exp", "disk", "plot", "subplot", "ntrt")) # ##dim 52192   13### 


##### get counts by plot year to crosscheck with excel file "Datasubset_CC Convergence###


da_origsiteyear<-da_orig %>%  distinct(field, year, exp, disk, plot, subplot, ntrt, .keep_all = T)

plotcountcheck<-da_origsiteyear%>%group_by(exp, field, year) %>%  tally()

plotcountcheck2<-da_origsiteyear%>%group_by(exp, year) %>%  tally()


plotcountcheck19822004<-subset(plotcountcheck, year < 2005)

sum(plotcountcheck19822004$n) #  6102 total plots-years 1982 to 2004 #### matches exactly w/ excel file## :-) 

sum(plotcountcheck$n)  ## 6396 total plot  years to be analyzed with full times series (note that this includes the E/W subplots that are compiled later...#### 



####Taxonomic and other small data fixes####

# Capitalize species to get rid of capitalization differences in spelling
da_orig$species <- toupper(as.character(da_orig$species))


###
da_orig$live <- 1
da_orig$sorted <- 1
da_orig$wood <- 0
da_orig$vasc <- 1

#Substitute outdated names with new names

da_orig$species <- gsub("APOCYNUM CANNABINUM", "APOCYNUM ANDROSAEMIFOLIUM", da_orig$species) ## MD 11/1 based off email with Eric##
da_orig$species <- gsub("MISC. FORB", "MISCELLANEOUS FORB", da_orig$species)
da_orig$species <- gsub("SEDGES", "CAREX SP.", da_orig$species)
da_orig$species <- gsub("QUERCUS RUBRUM", "QUERCUS RUBRA", da_orig$species)

# Set litter as not alive
sel<-da_orig$species == 'MISCELLANEOUS LITTER'
da_orig$live[sel] <- 0

sel<-grep("PINE", da_orig$species)
da_orig$live[sel] <- 0

# Code unsorted material as not being sorted
sel<-grep("MISCELLANEOUS", da_orig$species)
da_orig$sorted[sel] <- 0

sel<-grep("FUNGI", da_orig$species)
da_orig$sorted[sel] <- 0

sel<-grep("MOSS", da_orig$species)
da_orig$sorted[sel] <- 0

sel<-grep("LICHEN", da_orig$species)
da_orig$sorted[sel] <- 0

# Woody species
sel<-grep("ACER NEGUNDO", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("CEANOTHUS", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("CORYLUS", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("PINUS", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("PINE", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("POPULUS", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("QUERCUS", da_orig$species)
da_orig$wood[sel] <- 1
sel<-grep("ULMUS", da_orig$species)
da_orig$wood[sel] <- 1


##mistake w/ biomass for Asclepias syriaca##
da_orig["mass.above"][da_orig["mass.above"] == 1711.970] <- 1711.970/10


# Read in species attribute look up table
sp.df <- read.csv(here::here("data/Cedar_Creek_Plant_Taxon_List.csv"),header=T, fileEncoding = "UTF-8-BOM")



names(sp.df)[] <- tolower(names(sp.df))

#sp.df<-plyr::rename(sp.df, c("ï..species"="species"))
sp.df$species <- toupper(sp.df$species) 
sp.df$functional.group <- toupper(sp.df$functional.group)
sp.df$duration <- toupper(sp.df$duration)
sp.df$lifeform <- toupper(sp.df$lifeform)
sp.df$pathway <- toupper(sp.df$pathway)
sp.df$origin <- toupper(sp.df$origin)
sp.df$family <- toupper(sp.df$family)
names(sp.df[1])<-"species"



da_full <- merge(da_orig,sp.df, by='species', all.x=TRUE)

# Set Unknown species to missing
da_full$origin[da_full$origin == 'UNKNOWN'] <- NA
da_full$origin[da_full$origin == 'NATIVE AND/OR INTRODUCED'] <- NA


da_full$functional.group[da_full$functional.group == 'UNKNOW'] <- NA

da_full$duration[da_full$duration == 'UNKNOWN'] <- NA
da_full$duration[da_full$duration == 'BIENNIAL, PERENNIAL'] <- "BIENNIAL"

# Lump biennials in with annuals as they are pretty similar
# mostly weedy forbs 
unique(da_full$species[da_full$duration=="BIENNIAL"])
unique(da_full$species[da_full$duration=="BIENNIAL, PERENNIAL"])
unique(da_full$species[grep('BIENNIAL', da_full$duration)]) 

da_full$duration[grep('BIENNIAL', da_full$duration)] <- "ANNUAL"
unique(da_full$species[da_full$duration=="ANNUAL" & da_full$functional.group=="C4"])
unique(da_full$species[da_full$duration=="ANNUAL" & da_full$functional.group=="F"])

# Add in some species atrributes
sel <- da_full$species == "QUERCUS RUBRA"
da_full$functional.group[sel] <- "W"
da_full$lifeform[sel] <- "WOODY"
da_full$origin[sel] <- "NATIVE"

# Add in some species atrributes
sel <- da_full$species == "RUDBEKIA SEROTINA"
da_full$functional.group[sel] <- "F"
da_full$lifeform[sel] <- "FORB"
da_full$origin[sel] <- "NATIVE"

# Add in some species atrributes
sel <- da_full$species == "POA PRATENSIS"
da_full$origin[sel] <- "INTRODUCED"

sel <- da_full$species == "MISCELLANEOUS CAREX SP."
da_full$functional.group[sel] <- "S"
da_full$lifeform[sel] <- "SEDGE"
da_full$origin[sel] <- "NATIVE"


# Capitalize species to get rid of capitalization differences in spelling
da_full$species <- toupper(as.character(da_full$species))

##deal with some entries where sp. were weighed twice##
summary(freq <- ddply(da_full[da_full$live==1 & da_full$sorted==1, ], .(year, field, exp, plot, subplot, disk, ntrt, species), colwise(length, .(mass.above)))) ##takes a while##

freq$freq <- freq$mass.above

freq$mass.above <- NULL
freq[freq$freq > 1,]

doubles.df <- merge(freq[freq$freq > 1,], da_full[c("field", "exp","plot", "year", "species", "mass.above")], by=c("field", "exp","plot", "year", "species"))
doubles.df[c("field", "exp","plot", "year", "species", "mass.above")]

#####
# There are a few cases where there are multiple species weighed per sample. 
# Options are taking max, min, mean, or sum. 


# subset to live, sorted, herbaceous plants
d2 <- da_full[da_full$sorted ==1 & da_full$live ==1 & da_full$wood==0, c("field", "exp","plot", "subplot", "year", "disk", "ntrt",  "species", "mass.above")]

#subset that includes woody stuff####
d2woody <- da_full[da_full$sorted ==1 & da_full$live ==1, c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "species", "mass.above")]


#Add some columns for herbaceous + woody version ####
d3woody<-d2woody %>% mutate(expfieldplot=paste(exp,field,plot, sep = '_') , expntrtfieldyear= paste(exp,field, ntrt, year, sep = '_'), expntrtyear= paste(exp,ntrt, year, sep = '_'), ntrt2=ntrt) 
d3 <- d2 %>% mutate(expfieldplot=paste(exp,field,plot, sep = '_') , expntrtfieldyear= paste(exp,field, ntrt, year, sep = '_'), expntrtyear= paste(exp,ntrt, year, sep = '_'), ntrt2=ntrt) 


# remove woody species, species repeats, and seedling entries using new data file###
cleaned.d3woody <- read.csv(here::here("data/NEW.sp.decisions.csv"),header=T, fileEncoding = "UTF-8-BOM")
cleaned.d3woody <- cleaned.d3woody %>%
  filter(keep != "NO")

# combine sp. biomass with single species biomass 
species.keep <- cleaned.d3woody$species
cleaned.d3.sp <- d3woody %>%
  filter(species %in% species.keep)%>%
  unite("ID", 1:7, sep="_", remove = FALSE)

cleaned.d3 <-  left_join(cleaned.d3.sp, cleaned.d3woody)%>%
  mutate(species = recode(species, `APOCYNUM SP.` = "APOCYNUM ANDROSAEMIFOLIUM",
                          `ARISTIDA SP.` = "ARISTIDA BASIRAMEA",
                          `BROMUS SP.` = "BROMUS INERMIS",
                          `EQUISETUM SP.` = "EQUISETUM LAEVIGATUM",
                          `MELILOTUS SP.` = "MELILOTUS ALBA",
                          `OENOTHERA SP.` = "OENOTHERA BIENNIS",
                          `SILENE SP.` = "SILENE ANTIRRHINA"))

combine.mass <- cleaned.d3 %>% 
  dplyr::filter(keep == "COMBINE") %>%
  dplyr::group_by(ID, year, species)%>%
  dplyr::summarise(mass.above = sum(mass.above))%>%
  separate(ID, into = c("field", "exp", "plot", "subplot", "year", "disk", "ntrt"), remove=FALSE)
combine.mass$exp <- as.numeric(combine.mass$exp)
combine.mass$plot <- as.numeric(combine.mass$plot)
combine.mass$year <- as.numeric(combine.mass$year)
combine.mass$disk <- as.numeric(combine.mass$disk)
combine.mass$ntrt <- as.numeric(combine.mass$ntrt)

cleaned.d3 <- cleaned.d3 %>%
  filter(keep == "YES")%>%
  select(1:10)

cleaned.d3 <- rbind(cleaned.d3, combine.mass)
cleaned.d3 <- left_join(cleaned.d3, d3woody)%>%
  select(2:14)

# Transpose data to be a site by species matrix
da.wide.woody <- reshape(d3woody,
                        v.names="mass.above",
                        idvar=c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "expfieldplot", "expntrtyear"),
                        timevar="species",
                        direction="wide") 
da.wide.sp <- reshape(cleaned.d3.sp,
                         v.names="mass.above",
                         idvar=c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "expfieldplot", "expntrtyear"),
                         timevar="species",
                         direction="wide")
da.wide<- reshape(cleaned.d3,
                         v.names="mass.above",
                         idvar=c("field", "exp","plot", "subplot", "year", "disk", "ntrt", "expfieldplot", "expntrtyear"),
                         timevar="species",
                         direction="wide") 
# Fill in NA's with zeros in wide data set
da.wide[is.na(da.wide)] <- 0
da.wide.woody[is.na(da.wide.woody)] <- 0

##make nutrient and years factors#
da.wide$ntrt<-as.factor(da.wide$ntrt)
da.wide$year<-as.factor(da.wide$year)

da.wide.woody$ntrt<-as.factor(da.wide.woody$ntrt)
da.wide.woody$year<-as.factor(da.wide.woody$year)
##Combine east / west subplots into "whole" to match up with rest of data (from 2015)

# Set experiment and disturbance treatments as factors
da.wide$exp<-as.factor(da.wide$exp)
da.wide$disk<-as.factor(da.wide$disk)

da.wide.woody$exp<-as.factor(da.wide.woody$exp)
da.wide.woody$disk<-as.factor(da.wide.woody$disk)

#subset out the east/ west subplots plots and the whole plots
da.wideeastwest<-subset(da.wide, subplot=="East"|subplot=="West", row.names=NULL)
da.widewhole<-subset(da.wide, subplot=="Whole", row.names=NULL)

da.wideeastwest.woody<-subset(da.wide.woody, subplot=="East"|subplot=="West", row.names=NULL)
da.widewhole.woody<-subset(da.wide.woody, subplot=="Whole", row.names=NULL)

##Add the values from the subplots together##
EWaddall<-as.data.frame(da.wideeastwest%>%group_by(field,exp,plot,year,disk,ntrt,expfieldplot,expntrtyear, expntrtfieldyear, ntrt2) %>% summarise_if(is.numeric,mean)%>% mutate(subplot="Whole"))
EWaddall.woody<-as.data.frame(da.wideeastwest.woody%>%group_by(field,exp,plot,year,disk,ntrt,expfieldplot,expntrtyear, expntrtfieldyear, ntrt2) %>% summarise_if(is.numeric,mean)%>% mutate(subplot="Whole"))


##Merge the whole plots back together###
da.wide_allwhole<-rbind(da.widewhole,EWaddall)
da.wide_allwhole.woody<-rbind(da.widewhole.woody,EWaddall.woody)


dim(da.wide_allwhole) ## dimensions of the dataset = 6396  199##
dim(da.wide_allwhole.woody)
##reorder columns so that all the rows are in the same order across fields etc. 
da.wide5<-da.wide_allwhole %>% 
  arrange(year) %>% 
  arrange(plot) %>%
  arrange(exp) %>%
  arrange(field) 

da.wide5.woody<-da.wide_allwhole.woody %>% 
  arrange(year) %>% 
  arrange(plot) %>%
  arrange(exp) %>%
  arrange(field) 

### take out all plot years after 2004###

da.wide6<-subset(da.wide5, as.numeric(as.character(year)) < 2005)  ## dim 6102 rows ####
da.wide6.woody<-subset(da.wide5.woody, as.numeric(as.character(year)) < 2005)  ## dim 6102 rows ####

### tally up the years by each plot##
plotyearcheck<-da.wide6%>%group_by(exp, field, plot) %>%  tally()
plotyearcheck.woody<-da.wide6.woody%>%group_by(exp, field, plot) %>%  tally()
### for synchrony we only want plots that have at least 19 years of data (plots with original treatments the whole time series)

plotyearchecksynch<-plotyearcheck %>% mutate(expfieldplot=paste(exp,field,plot, sep = '_')) %>% mutate(includesynch=ifelse(n > 18, 1, 0))
plotyearchecksynch.woody<-plotyearcheck.woody %>% mutate(expfieldplot=paste(exp,field,plot, sep = '_')) %>% mutate(includesynch=ifelse(n > 18, 1, 0))
##dim 326 X 6 ####
plotyearchecksynchonly<-subset(plotyearchecksynch, includesynch==1)
plotyearchecksynchonly.woody<-subset(plotyearchecksynch.woody, includesynch==1)
#### dim 243 X 6 ###


da.widesynch= da.wide6 %>%  filter(expfieldplot %in% plotyearchecksynchonly$expfieldplot) 
da.widesynch.woody= da.wide6.woody %>%  filter(expfieldplot %in% plotyearchecksynchonly$expfieldplot) 

# dim = 5292 X 199 #####


# Return counts by plot and year to crosscheck with excel file "Datasubset_CC synchrony"


plotcountchecksych<-da.widesynch%>%group_by(exp, field, year) %>%  tally()
plotcountchecksych.woody<-da.widesynch.woody%>%group_by(exp, field, year) %>%  tally()

plotcountchecksynch2<-da.widesynch%>%group_by(exp, year) %>%  tally()  ## matches exactly w/ excel file##
plotcountchecksynch2.woody<-da.widesynch.woody%>%group_by(exp, year) %>%  tally()