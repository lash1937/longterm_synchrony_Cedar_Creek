#####################################
# This script is used to compare the component 
# metrics of synchrony and stability
########################

# Read in data and functions from source code
source(here::here("data_cleaning/subsetting_CC.R"))

#### Fig 1B. comparing population variance to community variance. ####
# prepare E001 data for analyses
sub_exp1_long <- da.widesynch %>% 
  dplyr::filter(exp == "1") %>% 
  tidyr::unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = FALSE) %>% 
  dplyr::select(-subplot, -ntrt2, -expntrtyear, -expntrtfieldyear)

sub_exp1_long <- sub_exp1_long %>% 
  tidyr::pivot_longer(9:length(sub_exp1_long),
                      names_to = "Species", values_to = "Abundance") %>% 
  filter(Abundance > 0)

sub_exp1_wide <- tidyr::pivot_wider(sub_exp1_long, names_from = year, 
                                    values_from = Abundance) %>% 
  select(-Species,-field,-exp,-plot,-disk,-ntrt,-expfieldplot)

sub_exp1_wide[is.na(sub_exp1_wide)] <- 0
plotnames<-c(unique(sub_exp1_long$uniqueID))
plotnames2<-str_sort(plotnames)

vreq_exp_1<-sub_exp1_wide

# create a blank matrix to fill 
matrix3<-matrix(NA, ncol=4, nrow=length(plotnames))
colnames(matrix3)<-c("plot", "classicVR", "comm", "pop")
plotnames_exp1<-c(unique(vreq_exp_1$uniqueID))

# loop through each unique ID (plot) and calculate pop variability,
# comm variability, and total variance ratio
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

# prepare E002 data for analyses
sub_exp2_long <- da.widesynch %>% 
  dplyr::filter(exp == "2") %>% 
  tidyr::unite(uniqueID, field, exp, plot, disk, ntrt, sep="_", remove = FALSE)%>% 
  dplyr::select(-subplot, -ntrt2, -expntrtyear, -expntrtfieldyear)

sub_exp2_long <- sub_exp2_long %>% tidyr::pivot_longer(9:length(sub_exp2_long),
                                                       names_to = "Species", 
                                                       values_to = "Abundance") %>% 
  filter(Abundance > 0)

sub_exp2_wide <- tidyr::pivot_wider(sub_exp2_long, names_from = year, 
                                    values_from = Abundance) %>% 
  select(-Species,-field,-exp,-plot,-disk,-ntrt,-expfieldplot)

sub_exp2_wide[is.na(sub_exp2_wide)] <- 0

plotnames<-c(unique(sub_exp2_long$uniqueID))
plotnames2<-str_sort(plotnames)

vreq_exp_2<-sub_exp2_wide
# create blank matrix to fill
matrix4<-matrix(NA, ncol=4, nrow=length(plotnames))  
colnames(matrix4)<-c("plot", "classicVR", "comm", "pop") 
plotnames_exp2<-c(unique(vreq_exp_2$uniqueID))

# Loop through each unique ID (plot) and calculate pop variability,
# comm variability, and total variance ratio
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

# combine E001 and E002 dataframes
popcomvar_exp12<-rbind(popcomvar_exp_1, popcomvar_exp_2)
popcomvar_exp12_2<- tidyr::separate(popcomvar_exp12, "plot", c("field", "exp", "plot", "disk", "ntrt"), sep="_")

# mutate nutrient column to reflect nitrogen concentration added
popcomvar_exp12_2<-popcomvar_exp12_2 %>% dplyr::mutate(Nitrogen=ntrt)
popcomvar_exp12_2$Nitrogen<-mapvalues(popcomvar_exp12_2$Nitrogen, 
                                      from=c( "1", "2", "3" ,"4", 
                                              "5", "6", "7", "8", "9"),
                                      to=c("0.0", "1.0", "2.0" ,"3.4",
                                           "5.4", "9.5", "17", "27.2", "0.0"))
popcomvar_exp12_2$Nitrogen<-as.numeric(as.character(popcomvar_exp12_2$Nitrogen))


# average across the disk and nutrients
popcomvar_exp12_2$disk<-as.factor(popcomvar_exp12_2$disk)
popcomvar_exp12_2$Nitrogen<-as.factor(popcomvar_exp12_2$Nitrogen)
popcomvar_exp12_2$classicVR<-as.numeric(popcomvar_exp12_2$classicVR)
popcomvar_exp12_2$comm<-as.numeric(popcomvar_exp12_2$comm)
popcomvar_exp12_2$pop<-as.numeric(popcomvar_exp12_2$pop)

# summarize mean population variance, and mean community variance
avgpopcommvar<- popcomvar_exp12_2 %>% dplyr::group_by(disk,Nitrogen) %>% 
  dplyr::summarize(meanVR = mean(classicVR), meanpop=mean(pop), meancomm=mean(comm))


my_virdis_pal <- c(viridis::viridis(n = 8, direction = -1))
my_virdis_pal[7:8] <- c("grey40", "black")

Fig1B<- ggplot()+
  geom_point(data = popcomvar_exp12_2, mapping =aes(x=pop,y=comm, col=Nitrogen,shape=disk),alpha=0.4)+
  geom_point(data= avgpopcommvar, mapping = aes(x=meanpop,y=meancomm, fill=Nitrogen,shape=disk),size=3)+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(0,0.6))+
  scale_y_continuous(limits=c(0,0.6))+
  labs(x="Aggregate Population \n Variability", y="Community Variability",tag = "B")+
  #coord_fixed()+
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 12, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = 'bottom')+
  geom_line(data = avgpopcommvar, mapping = aes(x=meanpop, y=meancomm, group=Nitrogen, col=Nitrogen))+
  annotate("text", x = 0.08, y=0.6, label = "Synchrony", color = "darkgrey", size = 3.5) + 
  annotate("text", x = 0.5, y=0.01, label = "Compensation", color = "darkgrey", size = 3.5) +  
  scale_shape_manual(name = "Disturbance",
                     labels = c("Intact", "Disturbed"),
                     values=c(21,24))+
  scale_fill_manual(#option = "D",direction=-1,
    na.value="grey72",
                       name="Nitrogen Addition",
                       breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                       labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"),
    values = my_virdis_pal)+
  scale_colour_manual(#option = "D",direction=-1,
    na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"),
    values = my_virdis_pal)+
  guides(fill=guide_legend(override.aes=list(shape=21), title = "Nitrogen Addition", title.position = "top", direction = "horizontal"))+ 
  guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "top", direction = "verticle"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"))) 

#### Fig 2B. comparing mean to the standard deviation of total biomass ####
# create total biomass data fram
biomass_df <- unique_ID_long %>%
  dplyr::group_by(uniqueID, year)%>%
  dplyr::summarize(total_biomass = sum(Abundance))

# find mean of total biomass and std dev of total biomass per plot
biomass_overtime <- biomass_df %>%
  dplyr::group_by(uniqueID)%>%
  dplyr::summarise(mean_biomass = mean(total_biomass), stdev_biomass = sd(total_biomass))%>%
  dplyr::mutate(hand_stab = mean_biomass/stdev_biomass)

biomass_all_sep <- tidyr::separate(biomass_overtime, "uniqueID", 
                                   c("field", "exp", "plot", "disk", "ntrt"), 
                                   sep="_", remove = FALSE)

biomass_all_cont<-biomass_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
biomass_all_cont$Nitrogen<-mapvalues(biomass_all_cont$Nitrogen, 
                                     from=c( "1", "2", "3" ,"4", 
                                             "5", "6", "7", "8", "9"),
                                     to=c("0.0", "1.0", "2.0" ,"3.4", 
                                          "5.4", "9.5", "17", "27.2", "0.0"))
biomass_all_cont$Nitrogen<-as.numeric(as.character(biomass_all_cont$Nitrogen))

#remove control treatment 9 -  no micronutrients
biomass_all_cont_minus9 <- biomass_all_cont%>%  
  dplyr::filter(ntrt < 9)

# remove outliers for plotting purposes
biomass_all_cont_minus9andoutliers <- biomass_all_cont_minus9

# calculate stability
st_all <- codyn::community_stability(unique_ID_long, 
                                     time.var = "year", 
                                     abundance.var = "Abundance", 
                                     replicate.var = "uniqueID")

# separate the unique ID string into its five identifiers
st_all_sep <- tidyr::separate(st_all, "uniqueID", 
                              c("field", "exp", "plot", "disk", "ntrt"), 
                              sep="_", remove = FALSE)

# mutate nitrogen treatments to reflect the concentrations applied
st_all_cont<-st_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
st_all_cont$Nitrogen<-mapvalues(st_all_cont$Nitrogen, 
                                from=c( "1", "2", "3" ,"4", 
                                        "5", "6", "7", "8", "9"),
                                to=c("0.0", "1.0", "2.0" ,"3.4", 
                                     "5.4", "9.5", "17", "27.2", "0.0"))
st_all_cont$Nitrogen<-as.numeric(as.character(st_all_cont$Nitrogen))

# subset out nutrient treatment 9, control plots without micronutrients
st_all_cont_minus9 <- st_all_cont%>%
  dplyr::filter(ntrt < 9)

# unite unique ID identifers in new dataset
st_all_cont_minus9 <- st_all_cont_minus9 %>% 
  tidyr::unite("uniqueID", 1:5, sep="_",remove = FALSE)

# soil disturbance treatment, disk, as a factor
st_all_cont_minus9$disk <- as.factor(st_all_cont_minus9$disk)

# combine dataframes
test_stab <- dplyr::left_join(biomass_all_cont_minus9, st_all_cont_minus9)

avg_test_stab <- test_stab %>%
  dplyr::group_by(disk, Nitrogen)%>%
  dplyr::summarise(mean_st = mean(stability))

# still calculate averages with outliers included
biomass_disk_N <- biomass_all_cont_minus9 %>%
  dplyr::group_by(disk, Nitrogen)%>%
  dplyr::summarise(meanofmean_biomass = mean(mean_biomass), meanofsd_biomass = mean(stdev_biomass))%>%
  dplyr::mutate(hand_stab = meanofmean_biomass/meanofsd_biomass)

# prepare data to plot
avg_test <- dplyr::left_join(avg_test_stab, biomass_disk_N)
biomass_disk_N$disk<-as.factor(biomass_disk_N$disk)
biomass_disk_N$Nitrogen<-as.factor(biomass_disk_N$Nitrogen)
biomass_all_cont_minus9andoutliers$disk<-as.factor(biomass_all_cont_minus9andoutliers$disk)
biomass_all_cont_minus9andoutliers$Nitrogen<-as.factor(biomass_all_cont_minus9andoutliers$Nitrogen)

# calculate a reference line for stability
{
  ref_stability <- 161.8789/62.76717 #control values
  x <- seq(0,500)
  y <- x*ref_stability
  ref_line <- as.data.frame(cbind(x,y))
  
  ggplot() +
    geom_line(data=ref_line, mapping=aes(x=x, y=y))
  
  Fig2B<- ggplot()+
    geom_point(data = biomass_all_cont_minus9andoutliers, mapping = 
                 aes(x=stdev_biomass,y=mean_biomass, 
                     col=Nitrogen,shape=disk),alpha=0.4)+
    geom_point(data = biomass_disk_N, mapping =
                 aes(x=meanofsd_biomass,y=meanofmean_biomass, 
                     fill=Nitrogen,shape=disk),size=3)+
    geom_abline(slope = ref_stability)+
    labs(x="Standard Deviation \n of Total Biomass", y="Mean of Total Biomass",tag = "B")+
    theme_bw() +
    #coord_fixed()+
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 14, angle = 90, hjust = .5, face = "plain"),
          legend.title = element_text(color = "black", size = 12, angle = 0, hjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12,angle = 0, hjust = 0, face = "plain"),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          legend.position = 'bottom')+
    geom_line(data = biomass_disk_N, mapping = 
                aes(x=meanofsd_biomass, y=meanofmean_biomass, 
                    group=Nitrogen, col=Nitrogen))+
    scale_shape_manual(name = "Disturbance",
                       labels = c("Intact", "Disturbed"),
                       values=c(21,24))+
    scale_fill_manual(#option = "D",direction=-1,
      na.value="grey72",
      name="Nitrogen Addition",
      breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
      labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"),
      values = my_virdis_pal)+
    scale_colour_manual(#option = "D",direction=-1,
      na.value="grey72",
      name="Nitrogen Addition",
      breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
      labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"),
      values = my_virdis_pal)+
    guides(fill=guide_legend(override.aes=list(shape=21), title = "Nitrogen Addition", title.position = "top", direction = "horizontal"))+ 
    guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "top", direction = "verticle"))+ 
    guides(theme(legend.title = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
                 legend.text=element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"))) 
  
  
}

# load ggplot objects from another script to create figure 1 AB and figure 2 AB
syncvsN <- readRDS(file = here::here("data/syncvsN.rds"))
stabvsN <- readRDS(file = here::here("data/stabvsN.rds"))

new_figure_1 <- syncvsN + Fig1B + plot_layout(ncol = 2, widths = c(2,1)) &
  theme(legend.position = 'none')

legend1A <- get_legend(syncvsN)
legend1B <- get_legend(Fig1B)
legends1 <- as_ggplot(legend1A) +  as_ggplot(legend1B) + plot_layout(ncol = 2, widths = c(1,2))

pdf(file="Figures/New_Figure1AB.pdf", width = 10, height = 6)
new_figure_1 / legends1 
dev.off()

new_figure_2 <- stabvsN + Fig2B + plot_layout(ncol = 2, widths = c(2,1)) &
  theme(legend.position = 'none')

legend2A <- get_legend(stabvsN)
legend2B <- get_legend(Fig2B)
legends2 <- as_ggplot(legend2A) +  as_ggplot(legend2B) + plot_layout(ncol = 2, widths = c(1,2))

pdf(file="Figures/New_Figure2AB.pdf", width = 10, height = 6)
new_figure_2 / legends2 
dev.off()

# modelling the effect of nitrogen and disk on population variability and community variability
popcomvar_exp12_2$Nitrogen <- as.numeric(as.character(popcomvar_exp12_2$Nitrogen))
# create unique grid variable
popcomvar_exp12_2 <- popcomvar_exp12_2 %>%
  mutate(
    grid = factor(paste0(field, exp)))

# log nitrogen for linear fit
popcomvar_exp12_2 <- popcomvar_exp12_2 %>%
  mutate(log_N = log(Nitrogen+1))

# test model with and without an interaction
pop_mod_int <- nlme::lme(pop ~ log_N*disk + field, random = (~1|grid), data=popcomvar_exp12_2)
pop_mod_noint <- nlme::lme(pop ~ log_N+disk + field, random = (~1|grid), data=popcomvar_exp12_2)
AIC(pop_mod_int, pop_mod_noint) # no interaction fit best

MuMIn::r.squaredGLMM(pop_mod_noint)
summary(pop_mod_noint)
anova(pop_mod_noint)

# test model with and without an interaction
comm_mod_int <- nlme::lme(comm ~ log_N*disk + field, random = (~1|grid), data=popcomvar_exp12_2)
comm_mod_noint <- nlme::lme(comm ~ log_N+disk + field, random = (~1|grid), data=popcomvar_exp12_2)
AIC(comm_mod_int,comm_mod_noint) # no interaction fit best
MuMIn::r.squaredGLMM(comm_mod_noint)
an.commmodnoint <- anova(comm_mod_noint)
summary(comm_mod_noint)
anova(comm_mod_noint)

# modelling the effect of nitrogen and disk on the mean and standard deviation of total biomass
# create unique grid variable
biomass_all_cont_minus9 <-biomass_all_cont_minus9 %>%
  mutate(
    grid = factor(paste0(field, exp)))
# log nitrogen for linear fit
biomass_all_cont_minus9 <- biomass_all_cont_minus9 %>%
  mutate(log_N = log(Nitrogen+1))

# test model with and without an interaction
mean_mod_int <- nlme::lme(mean_biomass~log_N*disk + field, random = (~1|grid), data = biomass_all_cont_minus9)
mean_mod_noint <- nlme::lme(mean_biomass~log_N + disk + field, random = (~1|grid), data = biomass_all_cont_minus9)
AIC(mean_mod_int,mean_mod_noint) # interaction fit best
MuMIn::r.squaredGLMM(mean_mod_int)
anova(mean_mod_int)
summary(mean_mod_int)
# test model with and without an interaction
std_mod_int <- nlme::lme(stdev_biomass~log_N*disk + field, random = (~1|grid), data = biomass_all_cont_minus9)
std_mod_noint <- nlme::lme(stdev_biomass~log_N+disk + field, random = (~1|grid), data = biomass_all_cont_minus9)
AIC(std_mod_int,std_mod_noint) # interaction fit best
MuMIn::r.squaredGLMM(std_mod_int)
anova(std_mod_int)
summary(std_mod_int)

# determine effect of N with and without disturbance, using N as a factor to get the effect at each N level
biomass_all_cont_factN <- biomass_all_cont_minus9 %>%
  mutate(fac_N = as.factor(biomass_all_cont_minus9$log_N))%>%
  mutate(fac_disk = as.factor(biomass_all_cont_minus9$disk))
std_mod_noint_fac <- nlme::lme(stdev_biomass ~  fac_N * fac_disk + field,
                          random = (~ 1 | grid), data = biomass_all_cont_factN)
EM_disk0_stdev <- emmeans::emmeans(std_mod_noint_fac, ~fac_N|fac_disk, at = list(disk = c(0)))#highN = log(27.2 +1)=3.3393220
pairs(EM_disk0_stdev)

mean_mod_noint_fac <- nlme::lme(mean_biomass ~  fac_N * fac_disk + field,
                               random = (~ 1 | grid), data = biomass_all_cont_factN)
EM_disk0_mean <- emmeans::emmeans(mean_mod_noint_fac, ~fac_N|fac_disk, at = list(disk = c(0)))#highN = log(27.2 +1)=3.3393220
pairs(EM_disk0_mean)

