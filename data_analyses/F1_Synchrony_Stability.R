
#####################################
# This script produces summary outputs of linear and quadratic models,
# which are tested to fit the relationships between global change drivers
# and synchrony and stability. The figures produced represent the effects 
# of nitrogen addition and soil disturbance on synchrony and stability.
#####################################

# Load libraries
library(ggeffects)
library(patchwork)
library(emmeans)

# Read in data and functions from source code
source(here::here("data_cleaning/subsetting_CC.R"))

#### Calculate Synchrony ####
# calculate the variance ratio for each plot, using 'year' as the time variable
VR_all <- codyn::variance_ratio(unique_ID_long, 
                                time.var = "year",
                                species.var = "Species", 
                                abundance.var = "Abundance",
                                bootnumber = 1, 
                                replicate = "uniqueID",
                                average.replicates = FALSE)

# separate unique ID string into it's five identifiers
VR_all_sep <- tidyr::separate(VR_all, "uniqueID", 
                              c("field", "exp", "plot", "disk", "ntrt"), 
                              sep="_")

# mutate nitrogen treatments to reflect the concentrations applied
VR_all_cont<-VR_all_sep %>% dplyr::mutate(Nitrogen=ntrt)
VR_all_cont$Nitrogen<-mapvalues(VR_all_cont$Nitrogen, 
                                from=c( "1", "2", "3" ,"4", 
                                        "5", "6", "7", "8", "9"),
                                to=c("0.0", "1.0", "2.0" ,"3.4", 
                                     "5.4", "9.5", "17", "27.2", "0.0"))
VR_all_cont$Nitrogen<-as.numeric(VR_all_cont$Nitrogen)

# subset out nutrient treatment 9, control plots without micronutrients
VR_all_cont_minus9 <- VR_all_cont%>% 
  dplyr::filter(ntrt < 9)

# unite unique ID identifers in new dataset
VR_all_cont_minus9 <- VR_all_cont_minus9 %>% 
  tidyr::unite("uniqueID", 1:5, sep="_",remove = FALSE)

# soil disturbance treatment, disk, as a factor
VR_all_cont_minus9$disk <- as.factor(VR_all_cont_minus9$disk)

#### Calculate Stability ####
# calculate stability for each plot using year as the time variable
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
st_all_cont$Nitrogen<-as.numeric(st_all_cont$Nitrogen)

# subset out nutrient treatment 9, control plots without micronutrients
st_all_cont_minus9 <- st_all_cont%>%
  dplyr::filter(ntrt < 9)

# unite unique ID identifers in new dataset
st_all_cont_minus9 <- st_all_cont_minus9 %>% 
  tidyr::unite("uniqueID", 1:5, sep="_",remove = FALSE)

# soil disturbance treatment, disk, as a factor
st_all_cont_minus9$disk <- as.factor(st_all_cont_minus9$disk)

#### Fig 1A. The effect of global change drivers on synchrony ####
# compare a linear model to a quadratic model to best describe the relationship 
# between synchrony and global change drivers

# create unique grid variable
VR_all_cont_minus9 <- VR_all_cont_minus9 %>%
  mutate(grid = factor(paste0(field, exp)))

# linear model 
mVRl_lme <- nlme::lme(VR ~  Nitrogen * disk + field,
                     random = (~ 1 | grid), data = VR_all_cont_minus9)
# quadratic model
mVRq_lme <- nlme::lme(VR ~  disk * poly(Nitrogen,2,raw=TRUE) + field, 
                      random = (~ 1 | grid), data= VR_all_cont_minus9)

# test model fit using AIC function
rawaic <- AIC(mVRl_lme, mVRq_lme)
nR <- dim(VR_all_cont_minus9)[1]
aictable(rawaic, nR) # linear model fit best

# determine the average trend across fields for plotting purposes
cfa_VR <- ggeffects::ggemmeans(mVRl_lme, terms=c("Nitrogen", "disk"))

# determine effect of disturbance with no N addition and high N addition
EM_controlN <- emmeans::emmeans(mVRl_lme, ~disk |Nitrogen, at = list(Nitrogen = c(0)))
pairs(EM_controlN)

EM_highN <- emmeans::emmeans(mVRl_lme, ~disk|Nitrogen, at = list(Nitrogen = c(27)))
pairs(EM_highN)

# plot new predicted lines to smooth the quadratic
Fig1A_newmod<- ggplot() +
  geom_point(data = VR_all_cont_minus9, aes(x=Nitrogen, y=VR, group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = cfa_VR, aes(x = x, y = predicted, 
                               group = group, color = group),linewidth = 1) +
  geom_ribbon(data = cfa_VR, aes(
    x = x,
    y = predicted,
    group= group,
    fill = group,
    ymin = conf.low,
    ymax = conf.high),
    alpha = 0.2,
    show.legend = F) +
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5,2.0))+
  scale_x_continuous(breaks = c(0, 10, 20, 30), labels=c("", "", "", ""))+
  ylab("Synchrony")+
  xlab("")+
  geom_hline(yintercept=1, color="darkgrey", linetype = "dashed") +
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, vjust = 0, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 14,
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) + labs(tag = "A")


# ---------------------------------------------------------------------------
#### Fig 1B. The effect of global change drivers on stability ####

# compare a linear model to a quadratic model to best describe the relationship 
# between stability and global change drivers

st_all_cont_minus9 <- st_all_cont_minus9 %>%
  mutate(grid = factor(paste0(field, exp)))

# linear model
mSTl_lme <-
  nlme::lme(stability ~ disk * Nitrogen + field, random = (~1 | grid), data= st_all_cont_minus9)
# quadratic model
mSTq_lme <-  nlme::lme(stability ~ disk * poly(Nitrogen,2,raw=TRUE) + field,random = (~1 | grid), 
                       data = st_all_cont_minus9)

# test model fit using AIC function
rawaic <- AIC(mSTl_lme, mSTq_lme)
nR <- dim(st_all_cont_minus9)[1]
aictable(rawaic, nR) # linear model fit best

# determine the average trend across fields for plotting purposes
cfa_ST <- ggeffects::ggemmeans(mSTl_lme, terms=c("Nitrogen", "disk"))

# plot new predicted lines
Fig1B_newmod<- ggplot() +
  geom_point(data = st_all_cont_minus9, aes(x=Nitrogen, y=stability, 
                                            group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = cfa_ST, aes(x = x, y = predicted, 
                                   group = group, color = group),linewidth = 1) +
  geom_ribbon(data = cfa_ST, aes(
    x = x,
    y = predicted,
    group= group,
    fill = group,
    ymin = conf.low,
    ymax = conf.high),
    alpha = 0.2,
    show.legend = F) +
  scale_colour_manual("legend", values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  ylab("Stability")+
  xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
  lims(y=c(.5,4.25))+
  labs(legend="Disturbance")+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, vjust = 0, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 14, 
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  theme(legend.position="top") + labs(tag = "B")

# produce final figure as a pdf
legend_fig1 <- get_legend(Fig1B_newmod)

quartz(width = 10, height = 5)

Figure1ab<- Fig1A_newmod / Fig1B_newmod + 
  plot_layout(ncol=1, widths = c(3,3)) & 
  theme(legend.position = "none")

pdf(file = "Figures/Figure1_AB.pdf",width = 4, height = 6)
Figure1ab / as_ggplot(legend_fig1) + plot_layout(heights=c(4,4,2))
dev.off()



