# Synchrony across nutrient treatment calculations

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
library(gifski)
library(png)
library(gganimate)
source(here("data_cleaning/cleaning_CC.R"))

calcSE <- function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

##starting file = exp12subset_2 dim(4092X 176)##

# want data formatted with species, year, location, abundance

unique_ID_exp12 <- unite(exp12subset_2, "uniqueID", c(field, exp, plot, disk, ntrt), sep="_")

#converting from wide to longform, keeps the zeros##

unique_ID_long <- pivot_longer(unique_ID_exp12, `mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`,
             names_to = "Species", values_to = "Abundance")

###just plot information#

st_all <- community_stability(unique_ID_long, time.var = "year", abundance.var = "Abundance", replicate.var = "uniqueID")

st_all_sep <- separate(st_all, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

st_all_cont<-st_all_sep %>% mutate(Nitrogen=ntrt)
st_all_cont$Nitrogen<-mapvalues(st_all_cont$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "9"),
                                to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "0.0"))
st_all_cont$Nitrogen<-as.numeric(as.character(st_all_cont$Nitrogen))

st_all_cont_minus9_first7 <- st_all_cont%>%
  filter(ntrt < 9)

##
st_all_cont_minus9 <- st_all_cont%>%
  filter(ntrt < 9)

lm_fit_st <- lm(stability ~ Nitrogen * disk, data=st_all_cont_minus9)
summary(lm_fit_st)
lm_fit_st <- lm(stability ~ Nitrogen +disk+Nitrogen*disk+field, data=st_all_cont_minus9)



summary(lm_fit_st)
anova(lm_fit_st)
plot(lm_fit_st)


hist(st_all_cont_minus9$stability)

st_all_cont_minus9$disk<-as.factor(st_all_cont_minus9$disk)


# graph length
p2 <- ggplot(data=st_all_cont_minus9, aes(x=Nitrogen, y=stability, group=disk, col=disk)) + 
  geom_jitter( size = 2,width = 0.2, height = 0, shape = 21)+
  geom_smooth(method = 'lm', size = 1.5)+
  #facet_wrap(.~field)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont_minus9, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="Nitrogen Addition", y="Community Stability") +
  #geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Undisturbed", "Disturbed"))


# let's calculate the variance ratio
# setting average.replicates to false allows us to
VR_all <- variance_ratio(unique_ID_long, time.var = "year",
                         species.var = "Species", abundance.var = "Abundance",
                         bootnumber = 1, replicate = "uniqueID",
                         average.replicates = FALSE)



VR_all_sep <- separate(VR_all, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

vr_m <- lm(VR~Nitrogen+disk+Nitrogen*disk+field, data = VR_all_cont_minus9)
summary(vr_m)
plot(vr_m)


anova(vr_m)
VR_averaged <- VR_all_cont %>% 
  group_by(disk) %>% 
  summarize(meanVR = mean(VR), seVR = calcSE(VR), CI=1.96*calcSE(VR))

# setting average.replicates to false allows us to

VR_all_cont<-VR_all_sep %>% mutate(Nitrogen=ntrt)
VR_all_cont$Nitrogen<-mapvalues(VR_all_cont$Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "9"),
                            to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "0.0"))
VR_all_cont$Nitrogen<-as.numeric(as.character(VR_all_cont$Nitrogen))

VR_all_cont_minus9 <- VR_all_cont%>%
  filter(ntrt < 9)

VR_averaged <- VR_all_cont_minus9 %>% 
  group_by(disk) %>% 
  summarize(meanVR = mean(VR), seVR = calcSE(VR), CI=1.96*calcSE(VR))

# let's try an example with just one of the unique ID's to make sure that we have 
# it right for a single treatment   

## plot specifications

#Create color palette
mycols <- colors()[c(578, 150)]
mypal <- palette(mycols)
names(mypal) = c("0", "1")
colScale <- scale_colour_manual(name = "Disturbance", values = mypal)


lm_fit <- lm(VR ~ Nitrogen * disk, data=VR_all_cont_minus9)
summary(lm_fit)


hist(VR_all_cont_minus9$VR)

VR_all_cont$disk<-as.factor(VR_all_cont$disk)
# graph length
p1 <- ggplot(data=VR_all_cont_minus9, aes(x=Nitrogen, y=VR, group=disk, col=disk)) +
  geom_jitter( size = 2,width = 0.2, height = 0, shape = 21)+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
                                      #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  #facet_wrap(.~field)+
  labs(x="Nitrogen Addition", y="Variance Ratio") +
  geom_hline(yintercept=1, color="darkgrey", linetype = "dashed") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5,2.0))+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                        name="Disturbance",
                        breaks=c("0", "1"),
                        labels=c("Undisturbed", "Disturbed"))




g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)



grid.arrange(g1, g2 , ncol = 1, heights = c(1, 1))



ggsave("vr_nut_disk.png", p1, dpi=300)

st_vr_all <- inner_join(st_all_cont_minus9, VR_all_cont_minus9)
p3 <- ggplot(data=st_vr_all, aes(x=VR, y=stability, group=disk, col=disk)) + 
  geom_jitter( size = 2,width = 0.2, height = 0, shape = 21)+
  geom_smooth(method = 'lm', size = 1.5)+
  ylim(0,5)+
  xlim(0,2)+
  geom_vline(xintercept = 1, color="grey")+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="Synchrony", y="Stability") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Undisturbed", "Disturbed"))

ggsave("stab_nut_disk.png", p2, dpi=300)

                                                                    
##### calculate stability --- fix with new subsetting method! ####

# stab_all_nozero <- community_stability(unique_ID_long_nozero, time.var = "year",
#                                      abundance.var = "Abundance", replicate.var = "uniqueID")
# stab_all_sep_nozero <- separate(stab_all_nozero, "uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_")
# 
# stab_averaged <- stab_all_sep_nozero %>% group_by(disk,ntrt) %>%
#   summarize(meanST = mean(stability), seST = calcSE(stability), CIST=1.96*calcSE(stability))
# 
# 
# stab_averaged$ntrt2 <- stab_averaged$ntrt
# stab_averaged$ntrt2  <- mapvalues(stab_averaged$ntrt2, from=c("9", "1", "2", "4",  "6"),
#                                        to=c("Control", "Control+micro", "Low N + micro", "Med N + micro", "High N + micro"))
# stab_averaged$ntrt2 <- factor(stab_averaged$ntrt2, levels = c("Control", "Control+micro", "Low N + micro", "Med N + micro", "High N + micro"))
# 
# ggplot() +
#   facet_wrap(~disk, ncol=1)+
#   geom_point(data=stab_averaged, aes(x=ntrt2, y=meanST)) +
#   #geom_line() +
#   geom_errorbar(data=stab_averaged, aes(x= ntrt2, ymin = meanST - CIST,
#                                          ymax = meanST + CIST, width=0.2)) +
#   labs(x="Nitrogen Treatment", y="Stability") +
#   geom_hline(yintercept=1, color="darkgrey")+
#   theme_classic()



## sub-setting timescale specific breakdown into exp 1 and exp 2 to handle missing data
# experiment 1
sub_exp1_long <- exp12subset_2 %>% filter(exp == "1") %>% unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = TRUE)%>% 
  select(-subplot,-nadd, -ntrt2, -repwithintrt, -expntrtyear)%>%
  pivot_longer(`mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`,
               names_to = "Species", values_to = "Abundance")%>%filter(Abundance > 0)

sub_exp1_wide <- pivot_wider(sub_exp1_long, names_from = year, values_from = Abundance)%>%select(-Species)

sub_exp1_wide[is.na(sub_exp1_wide)] <- 0

plotnames<-c(unique(sub_exp1_long$uniqueID))

plotnames2<-str_sort(plotnames)

matrix1<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix1)<-c("plot", "classic", "short", "long")

vreq_exp_1<-sub_exp1_wide

for(s in 1:length(plotnames2)){
  current.plot<-plotnames2[s]
  temp<-vreq_exp_1 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  timescale_s <-aggts(timescale,timescale$ts[timescale$ts<4])
  timescale_l <-aggts(timescale,timescale$ts[timescale$ts>=4])
  matrix1[s,2] <- vr[[3]]
  matrix1[s,1] <- current.plot
  matrix1[s,3] <- timescale_s[[3]]
  matrix1[s,4] <- timescale_l[[3]]
  
}

TSVR_exp_1 <- as.data.frame(matrix1)
TSVR_exp_1 <- TSVR_exp_1 %>% pivot_longer(cols = 2:4, names_to = "VR_type", values_to = "value")%>%
  separate("plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

TSVR_exp_1$value <- as.numeric(as.character(TSVR_exp_1$value))

averagenutTSVR_exp_1<- TSVR_exp_1 %>% group_by(disk,ntrt,VR_type) %>% 
  summarize(meanVR = mean(value), seVR = calcSE(value), CI=1.96*calcSE(value))

# experiment 2
sub_exp2_long <- exp12subset_2 %>% filter(exp == "2")%>% unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = TRUE)%>%
  select(-subplot,-nadd, -ntrt2, -repwithintrt, -expntrtyear)%>%
  pivot_longer(`mass.above.ACHILLEA MILLEFOLIUM(LANULOSA)`:`mass.above.VIOLA SP.`,
               names_to = "Species", values_to = "Abundance")%>%filter(Abundance > 0)

sub_exp2_wide <- pivot_wider(sub_exp2_long, names_from = year, values_from = Abundance) %>% select(-'2003', -Species)

sub_exp2_wide[is.na(sub_exp2_wide)] <- 0

sub_exp2_wide2<-sub_exp2_wide[rowSums(sub_exp2_wide[,2:20])>0,]

plotnames<-c(unique(sub_exp2_long$uniqueID))

plotnames2<-str_sort(plotnames)

matrix2<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix2)<-c("plot", "classic", "short", "long")


vreq_exp_2<-sub_exp2_wide2

for(s in 1:length(plotnames2)){
  current.plot<-plotnames2[s]
  temp<-vreq_exp_2 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  timescale_s <-aggts(timescale,timescale$ts[timescale$ts<4])
  timescale_l <-aggts(timescale,timescale$ts[timescale$ts>=4])
  matrix2[s,2] <- vr[[3]]
  matrix2[s,1] <- current.plot
  matrix2[s,3] <- timescale_s[[3]]
  matrix2[s,4] <- timescale_l[[3]]
  
}

TSVR_exp_2 <- as.data.frame(matrix2)
TSVR_exp_2 <- TSVR_exp_2 %>% pivot_longer(cols = 2:4, names_to = "VR_type", values_to = "value")%>%
  separate("plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

TSVR_exp_2$value <- as.numeric(as.character(TSVR_exp_2$value))

averagenutTSVR_exp_2<- TSVR_exp_2 %>% group_by(disk,ntrt,VR_type) %>% 
  summarize(meanVR = mean(value), seVR = calcSE(value), CI=1.96*calcSE(value))

TSVR_df <- rbind(averagenutTSVR_exp_1, averagenutTSVR_exp_2)

## visualize TSVR results ##


disturbance <- c("0" = "Undisturbance", "1" = "Disturbance")
ggplot() + 
  facet_wrap(~disk, ncol=1, labeller = labeller(disk = disturbance))+
  geom_point(data=TSVR_df, aes(x=ntrt, y=meanVR, group=VR_type, color=VR_type)) + 
  #geom_line() +
  geom_errorbar(data=TSVR_df, aes(x= ntrt, ymin = meanVR - CI, 
                                         ymax = meanVR + CI, group=VR_type, color=VR_type), width=0.2) +
  labs(x="Nitrogen Treatment", y="Variance Ratio") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_classic()+
  scale_colour_discrete(name="Time\nScale",
                    breaks=c("classic", "long", "short"),
                    labels=c("Classic", "Long", "Short"))
  


#### matrix for all aggregate properties exp 1 ####

vreq_exp_1<-sub_exp1_wide

# MCS find rows and or columns that are all 0s
empty <- c()
temp <- vreq_exp_1[,c(2:24)]
# looping the rows
for (i in 1:nrow(temp)){
  # counter for blank values in 
  # each row
  count = 0
  # looping through columns
  for(j in 1:ncol(temp)){
    # checking if the value is blank
    if(isTRUE(temp[i,j] == 0)){
      count = count + 1
    }
  }
  # if count is equivalent to number 
  # of columns
  if(count == ncol(temp)){
    # append row number
    empty <- append(empty,i)
  }
}# empty = NULL, all good.

matrix3<-matrix(NA, ncol=8, nrow=length(plotnames))
colnames(matrix3)<-c("plot", "classic", "cshort", "clong", "com_short", "com_long", "pop_short", "pop_long" )
plotnames_exp1<-c(unique(vreq_exp_1$uniqueID))

for(s in 1:length(plotnames_exp1)){
  current.plot<-plotnames_exp1[s]
  temp<-vreq_exp_1 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  timescale_s <-aggts(timescale,timescale$ts[timescale$ts<4])
  timescale_l <-aggts(timescale,timescale$ts[timescale$ts>=4])
  matrix3[s,2] <- vr[[3]]
  matrix3[s,1] <- current.plot
  matrix3[s,3] <- timescale_s[[3]]
  matrix3[s,4] <- timescale_l[[3]]
  matrix3[s,5] <- timescale_s[[1]]
  matrix3[s,6] <- timescale_l[[1]]
  matrix3[s,7] <- timescale_s[[2]]
  matrix3[s,8] <- timescale_l[[2]]
}




matrix3<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix3)<-c("plot", "classicVR", "comm", "pop")
plotnames_exp1<-c(unique(vreq_exp_1$uniqueID))


for(s in 1:length(plotnames_exp1)){
  current.plot<-plotnames_exp1[s]
  temp<-vreq_exp_1 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  matrix3[s,2] <- vr[[3]]
  matrix3[s,1] <- current.plot
  matrix3[s,3]<- vr[[1]]
  matrix3[s,4]<- vr[[2]]
}


popcomvar_exp_1 <- as.data.frame(matrix3)


#####
popcomevar_exp_1_wide <- TSVR_exp_3 %>% pivot_longer(cols = 2:4, names_to = "time", values_to = "value")%>%
  separate("plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

TSVR_exp_3$time <- as.factor(TSVR_exp_3$time)
TSVR_exp_3$VR_type <- NA
TSVR_exp_3$VR_type[TSVR_exp_3$time %in% c("classic","cshort","clong")] <- c("weighted")
TSVR_exp_3$VR_type[TSVR_exp_3$time %in% c("com_short","com_long")] <- c("CV_comm")
TSVR_exp_3$VR_type[TSVR_exp_3$time %in% c("pop_short","pop_long")] <- c("CV_ip")

TSVR_exp_3$timescale <- NA
TSVR_exp_3$timescale[TSVR_exp_3$time %in% c("classic")] <- c("classic")
TSVR_exp_3$timescale[TSVR_exp_3$time %in% c("cshort","com_short","pop_short")] <- c("short")
TSVR_exp_3$timescale[TSVR_exp_3$time %in% c("clong", "com_long","pop_long")] <- c("long")

TSVR_exp_3 <- TSVR_exp_3 %>% select(-time)
TSVR_exp_3 <- TSVR_exp_3 %>% unite(uniqueID, 1:5, sep= "_", remove = FALSE)

# want all nutrients
TSVR_exp_3_wide <- TSVR_exp_3 %>%
  pivot_wider(names_from = VR_type, values_from = value)%>%
  filter(timescale == c("short", "long"))

# control <- c("1" , "9")
# TSVR_exp_3_wide <- TSVR_exp_3 %>%
#   pivot_wider(names_from = VR_type, values_from = value)%>%
#   filter(timescale == c("short", "long"), ntrt %in% control)

## classic variance ratio for control plots combined ##

classic_control <- TSVR_exp_1 %>%
  filter(VR_type == "classic", ntrt %in% control)

classic_plot <- ggplot(classic_control, aes(x=VR_type, y=value))+geom_boxplot()+
  ylab("classic variance ratio")+
  xlab("")+
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 12))+
  geom_abline(intercept = 1, slope = 0, linetype = "dotted")



#### matrix for all aggregate properties exp 2 ####

matrix4<-matrix(NA, ncol=8, nrow=length(plotnames))
colnames(matrix4)<-c("plot", "classic", "cshort", "clong", "com_short", "com_long", "pop_short", "pop_long" )

vreq_exp_2<-sub_exp2_wide2
plotnames_exp2 <- c(unique(vreq_exp_2$uniqueID))

for(s in 1:length(plotnames_exp2)){
  current.plot<-plotnames_exp2[s]
  temp<-vreq_exp_2 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  timescale_s <-aggts(timescale,timescale$ts[timescale$ts<4])
  timescale_l <-aggts(timescale,timescale$ts[timescale$ts>=4])
  matrix4[s,2] <- vr[[3]]
  matrix4[s,1] <- current.plot
  matrix4[s,3] <- timescale_s[[3]]
  matrix4[s,4] <- timescale_l[[3]]
  matrix4[s,5] <- timescale_s[[1]]
  matrix4[s,6] <- timescale_l[[1]]
  matrix4[s,7] <- timescale_s[[2]]
  matrix4[s,8] <- timescale_l[[2]]
  
}

# make into dataframe 
TSVR_exp_4 <- as.data.frame(matrix4)
TSVR_exp_4 <- TSVR_exp_4 %>% pivot_longer(cols = 2:8, names_to = "time", values_to = "value")%>%
  separate("plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

TSVR_exp_4$time <- as.factor(TSVR_exp_4$time)
TSVR_exp_4$VR_type <- NA
TSVR_exp_4$VR_type[TSVR_exp_4$time %in% c("classic","cshort","clong")] <- c("weighted")
TSVR_exp_4$VR_type[TSVR_exp_4$time %in% c("com_short","com_long")] <- c("CV_comm")
TSVR_exp_4$VR_type[TSVR_exp_4$time %in% c("pop_short","pop_long")] <- c("CV_ip")

TSVR_exp_4$timescale <- NA
TSVR_exp_4$timescale[TSVR_exp_4$time %in% c("classic")] <- c("classic")
TSVR_exp_4$timescale[TSVR_exp_4$time %in% c("cshort","com_short","pop_short")] <- c("short")
TSVR_exp_4$timescale[TSVR_exp_4$time %in% c("clong", "com_long","pop_long")] <- c("long")

TSVR_exp_4 <- TSVR_exp_4 %>% select(-time)
TSVR_exp_4 <- TSVR_exp_4 %>% unite(uniqueID, 1:5, sep= "_", remove = FALSE)

# want all nutrients
TSVR_exp_4_wide <- TSVR_exp_4 %>%
  pivot_wider(names_from = VR_type, values_from = value)%>%
  filter(timescale == c("short", "long"))

# Dataframe for both experiments ####
TSVR_all <- rbind(TSVR_exp_3, TSVR_exp_4)
TSVR_all_wide <- TSVR_all %>%
  pivot_wider(names_from = VR_type, values_from = value)%>%
  filter(timescale == c("short", "long"))

# make it so that short comes first then long
TSVR_all_wide$timescale <- factor(TSVR_all_wide$timescale, levels = c("short","long"))

# name experiments so figs are information
TSVR_all_wide$exp <- factor(TSVR_all_wide$exp, labels = c('1' = "not disturbed", '2' = "disturbed"))

# ntrt levels from character to numeric so figure can have gradient fill
TSVR_all_wide$ntrt <- as.numeric(TSVR_all_wide$ntrt)

# create nutrient groups
TSVR_all_wide$ntrt_group <- NA
TSVR_all_wide$ntrt_group[TSVR_all_wide$ntrt %in%c(2,3)] <- 'low N'
TSVR_all_wide$ntrt_group[TSVR_all_wide$ntrt%in%c(4,5,6)] <- 'high N'
TSVR_all_wide$ntrt_group[TSVR_all_wide$ntrt%in%c(1,9)] <- 'control N'

# order groups
TSVR_all_wide$ntrt_group <- factor(TSVR_all_wide$ntrt_group, levels = c('control N', 'low N', 'high N'))

# community and population variability
matrix5<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix5)<-c("plot", "classic", "cvcomm", "cvcommip")

vcomm_exp_2<-sub_exp2_wide2
plotnames_exp2 <- c(unique(vcomm_exp_2$uniqueID))

for(s in 1:length(plotnames_exp2)){
  current.plot<-plotnames_exp2[s]
  temp<-vreq_exp_2 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  matrix5[s,2] <- vr[[3]]
  matrix5[s,3] <- vr[[1]]
  matrix5[s,4] <- vr[[2]]
  matrix5[s,1] <- current.plot
}

matrix6<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix6)<-c("plot", "classic", "cvcomm", "cvcommip")

plotnames_exp1<-c(unique(vreq_exp_1$uniqueID))

for(s in 1:length(plotnames_exp1)){
  current.plot<-plotnames_exp1[s]
  temp<-vreq_exp_1 %>% 
    filter(uniqueID %in% current.plot) %>% 
    select(-uniqueID)
  temp<-as.matrix(temp)
  vr<-vreq_classic(temp)
  timescale<-tsvreq_classic(temp)
  matrix6[s,2] <- vr[[3]]
  matrix6[s,3] <- vr[[1]]
  matrix6[s,4] <- vr[[2]]
  matrix6[s,1] <- current.plot
}


comm_v2 <- as.data.frame(matrix5)
comm_v1 <- as.data.frame(matrix6)
comm_var <- rbind(comm_v1, comm_v2)
comm_var <- comm_var %>%
  separate("plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")%>%
  mutate(Nitrogen = case_when( ntrt == 1 ~ "0.0",
                               ntrt == 2 ~ "1.0",
                               ntrt == 3 ~ "2.0",
                               ntrt == 4 ~ "3.4",
                               ntrt == 5 ~ "5.4",
                               ntrt == 6 ~ "9.5",
                               ntrt == 9 ~ "0.0"))
comm_var$Nitrogen<-as.numeric(as.character(comm_var$Nitrogen))
comm_var$cvcomm<-as.numeric(as.character(comm_var$cvcomm))
comm_var$cvcommip<-as.numeric(as.character(comm_var$cvcommip))
commvar <- ggplot(data=comm_var, aes(x=Nitrogen, y=cvcomm, group=disk, col=disk)) + 
  geom_point( size = 2, position='jitter')+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="Nitrogen Addition", y="Community Variability") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                      name="",
                      breaks=c("0", "1"),
                      labels=c("Not Disturbed", "Disturbed"))

popvar <- ggplot(data=comm_var, aes(x=Nitrogen, y=cvcommip, group=disk, col=disk)) + 
  geom_point( size = 2, position='jitter')+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="Nitrogen Addition", y="Aggregate Population Variability") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                      name="",
                      breaks=c("0", "1"),
                      labels=c("Not Disturbed", "Disturbed"))




lm_ipvar <- lm(cvcommip ~ Nitrogen * disk, data=comm_var)
summary(lm_ipvar)

lm_commvar <- lm(cvcomm ~ Nitrogen * disk, data=comm_var)
summary(lm_commvar)

sync_decomp <- comm_var %>%
  group_by(disk, Nitrogen)%>%
  summarize(avg_pop = mean(cvcommip), avg_comm = mean(cvcomm))

sync_decom_disk <- sync_decomp %>%
  group_by(Nitrogen)%>%
  summarize(avg_pop_diff = diff(avg_pop), avg_comm_diff = diff(avg_comm))

sync_decomp_N <- sync_decomp %>%
  group_by(disk)%>%
  summarize(avg_pop_diff = diff(avg_pop), avg_comm_diff = diff(avg_comm))


commvar_stab1 <- community_stability(unique_ID_long,  time.var = "year", abundance.var = "Abundance", replicate.var = "uniqueID" )
commvar_stab1 <- commvar_stab1%>% separate("uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_", remove = FALSE)%>%
  mutate(Nitrogen = case_when( ntrt == 1 ~ "0.0",
                               ntrt == 2 ~ "1.0",
                               ntrt == 3 ~ "2.0",
                               ntrt == 4 ~ "3.4",
                               ntrt == 5 ~ "5.4",
                               ntrt == 6 ~ "9.5",
                               ntrt == 9 ~ "0.0"))
commvar_stab1$Nitrogen<-as.numeric(as.character(commvar_stab1$Nitrogen))

stability_commvar <- inner_join(commvar_stab1,comm_var)

popstab <- ggplot(data=stability_commvar, aes(x=cvcommip, y=stability,  shape=disk, col=Nitrogen )) + 
  geom_point( size = 2, position='jitter')+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="aggregate population variability", y="stability") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())


commstab <- ggplot(data=stability_commvar, aes(x=cvcomm, y=stability,  col=Nitrogen, shape=disk)) + 
  geom_point( size = 2, position='jitter')+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="community variability", y="stability") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())


popcomm <- ggplot(data=comm_var, aes(x=cvcomm, y=cvcommip, shape=disk, col=Nitrogen)) + 
  geom_point( size = 2, position='jitter')+
  geom_smooth(method = 'lm', size = 1.5)+
  #geom_point(data=VR_all_cont, aes(x=Nitrogen, y=VR, group=disk, col=disk), size = 3) + 
  #geom_errorbar(data=VR_all_cont, aes(x= Nutrient, ymin = meanVR - CI, 
  #ymax = meanVR + CI, group=disk, color=disk), width=0.2) +
  labs(x="Community Variability", y="Aggregate Population Variability") +
  geom_hline(yintercept=1, color="darkgrey") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())





comm_var_0 <- comm_var %>%
  filter(disk == 0)
sync_decomp0 <- sync_decomp %>%
  filter(disk == 0)

sync_decomp$Nitrogen <- as.factor(sync_decomp$Nitrogen)
comm_var$Nitrogen <- as.factor(comm_var$Nitrogen)
comm_var_0$Nitrogen <- as.factor(comm_var_0$Nitrogen)
sync_decomp0$Nitrogen <- as.factor(sync_decomp0$Nitrogen)
sync_decomp$disk <- as.factor(sync_decomp$disk)
sync_decomp0$disk <- as.factor(sync_decomp0$disk)
comm_var$disk <- as.factor(comm_var$disk)
comm_var_0$disk <- as.factor(comm_var_0$disk)
p8<- ggplot()+
  geom_point(data = comm_var_0, mapping =aes(x=cvcommip,y=cvcomm, col=Nitrogen,shape=disk),alpha=0.4,size=0)+
  geom_point(data= sync_decomp0, mapping = aes(x=avgpop,y=avgcomm, fill=Nitrogen,shape=disk),size=0)+
  #transition_reveal(cvcommip)+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(0,.5))+
  scale_y_continuous(limits=c(0,.5))+
  labs(x="Aggregate Population Variability", y="Community Variability")+
  coord_fixed()+
  theme_classic2()+
  #geom_line(data = sync_decomp0, mapping = aes(x=avgpop, y=avgcomm, group=Nitrogen, col=Nitrogen))+
  scale_shape_manual(name = "Disturbance",
                     labels = c("Undisturbed", "Disturbed"),
                     values=c(21,24))+
  scale_fill_brewer(palette = "Paired",
                   name="Nitrogen Addition",
                  breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                  labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  scale_color_brewer(palette = "Paired",
                    name="Nitrogen Addition",
                    breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                     labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
ggsave("blank.png", p8, dpi=300)


p9<- ggplot()+
  geom_point(data = comm_var_0, mapping =aes(x=cvcommip,y=cvcomm, col=Nitrogen,shape=disk),alpha=0.4)+
  geom_point(data= sync_decomp0, mapping = aes(x=avgpop,y=avgcomm, fill=Nitrogen,shape=disk),size=3)+
  #transition_reveal(cvcommip)+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(0,.5))+
  scale_y_continuous(limits=c(0,.5))+
  labs(x="Aggregate Population Variability", y="Community Variability")+
  coord_fixed()+
  theme_classic2()+
  geom_line(data = sync_decomp0, mapping = aes(x=avgpop, y=avgcomm, group=Nitrogen, col=Nitrogen))+
  scale_shape_manual(name = "Disturbance",
                     labels = c("Undisturbed", "Disturbed"),
                     values=c(21,24))+
  scale_fill_brewer(palette = "Paired",
                    name="Nitrogen Addition",
                    breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                    labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  scale_color_brewer(palette = "Paired",
                     name="Nitrogen Addition",
                     breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                     labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  guides(fill=guide_legend(override.aes=list(shape=21))) 

ggsave("pop_comm_0.png", p9, dpi=300)


p10<- ggplot()+
  geom_point(data = comm_var, mapping =aes(x=cvcommip,y=cvcomm, col=Nitrogen,shape=disk),alpha=0.4)+
  geom_point(data= sync_decomp, mapping = aes(x=avg_pop,y=avg_comm, fill=Nitrogen,shape=disk),size=3)+
  #transition_reveal(cvcommip)+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(0,.5))+
  scale_y_continuous(limits=c(0,.5))+
  labs(x="Aggregate Population Variability", y="Community Variability")+
  coord_fixed()+
  theme_classic2()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  geom_line(data = sync_decomp, mapping = aes(x=avg_pop, y=avg_comm, group=Nitrogen, col=Nitrogen))+
  scale_shape_manual(name = "Disturbance",
                     labels = c("Undisturbed", "Disturbed"),
                     values=c(21,24))+
  scale_fill_brewer(palette = "Paired",
                    name="Nitrogen Addition",
                    breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                    labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  scale_color_brewer(palette = "Paired",
                     name="Nitrogen Addition",
                     breaks=c(0, 1, 2, 3.4, 5.4, 9.5),
                     labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+ guides(shape = guide_legend(override.aes = list(size = 0.5)))+theme(legend.title = element_text( size=10), legend.text=element_text(size=8))

ggsave("pop_com.png", p10, dpi=300)



# visualize synchrony and stability faceted by nutrient concentration
sync_stab_minus9 <- full_join(VR_all, st_all, by = "uniqueID")%>%
  select(-c("lowerCI", "upperCI", "nullmean"))%>%
  separate("uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_")%>%
  mutate(Nitrogen=ntrt)

sync_stab_minus9 $Nitrogen<-mapvalues(sync_stab_minus9 $Nitrogen, from=c( "1", "2", "3" ,"4", "5", "6", "9"),
                                 to=c("0.0", "1.0", "2.0" ,"3.4", "5.4", "9.5", "0.0"))
sync_stab_minus9 $Nitrogen<-as.numeric(as.character(sync_stab_minus9 $Nitrogen))
# remove control w/o micronutrients
sync_stab_minus9  <- sync_stab_minus9 %>%
  filter(ntrt < 9)

sync_stab_avg <- sync_stab_minus9 %>%
  group_by(disk,Nitrogen)%>%
  summarize(meanVR = mean(VR), meanST = mean(stability))

sync_stability <- ggplot(sync_stab_minus9, aes(x=VR, y=stability, group=disk, col=disk))+
  geom_point( size = 1)+
  geom_smooth(method = 'lm', size = 1.5)+
  labs(x="Synchrony", y="Stability") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 13, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 13, angle = 90, hjust = .5, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_colour_manual(values = c("#D55E00", "#0072B2"),
                      name="",
                      breaks=c("0", "1"),
                      labels=c("Not Disturbed", "Disturbed"))+
  facet_wrap(~ Nitrogen, ncol = 2)




  
