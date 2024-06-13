
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
# compare a linear model to a quadratic model (N and log(N)) to best describe the relationship 
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
VR_all_cont_minus9 <- VR_all_cont_minus9 %>%
  mutate(log_N = log(Nitrogen+1))

# log linear model 
mVRl_lme_log <- nlme::lme(VR ~  log_N * disk + field,
                      random = (~ 1 | grid), data = VR_all_cont_minus9)
# log quadratic model
mVRq_lme_log <- nlme::lme(VR ~  disk * poly(log_N,2,raw=TRUE) + field, 
                      random = (~ 1 | grid), data= VR_all_cont_minus9)

# test  model fits using AIC function
rawaic <- AIC(mVRl_lme, mVRq_lme, mVRl_lme_log, mVRq_lme_log)
nR <- dim(VR_all_cont_minus9)[1]
aictable(rawaic, nR) # log linear model fit best

MuMIn::r.squaredGLMM(mVRl_lme_log)
an.mVRl_log <- anova(mVRl_lme_log)

# determine the average trend across fields for plotting purposes
cfa_VR <- ggeffects::ggemmeans(mVRl_lme_log, terms=c("log_N", "disk")) %>% 
  as_tibble() %>%
  mutate(xo = expm1(x))

# determine effect of disturbance with no N addition and high N addition

em.mod <- emtrends(mVRl_lme_log, "disk", var = "log_N")

EM_controlN <- emmeans::emmeans(mVRl_lme_log, ~disk |log_N, at = list(log_N = c(0.0))) #control = log(0+1)=0.0
pairs(EM_controlN)

EM_highN <- emmeans::emmeans(mVRl_lme_log, ~disk|log_N, at = list(log_N = c(3.3393220)))#highN = log(27.2 +1)=3.3393220
pairs(EM_highN)
summary(mVRl_lme_log)
anova(mVRl_lme_log)
# plot new predicted lines to smooth the quadratic
Fig1A_newmod<- ggplot() +
  geom_point(data = VR_all_cont_minus9, aes(x=Nitrogen, y=VR, group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = cfa_VR, aes(x = xo, y = predicted, 
                               group = group, color = group),linewidth = 1) +
  geom_ribbon(data = cfa_VR, aes(
    x = xo,
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
                      labels=c("Intact", "Disturbed"))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5,2.0))+
  scale_x_continuous(breaks = c(0, 1, 2, 3.4, 5.4, 9.5, 17.0, 27.2), labels = c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2")) +
  ylab("Synchrony")+
  xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
  geom_hline(yintercept=1, color="darkgrey", linetype = "dashed") +
  theme_bw()+
  #coord_fixed()+
  annotate("text", x = 1.2, y=1.8, label = "Synchrony", color = "darkgrey", size = 3.5) + 
  annotate("text", x = 2.0, y=0.1, label = "Compensation", color = "darkgrey", size = 3.5) +  
  theme(axis.text.x = element_text(color = "grey20", size = 12,
                                   angle = 45, hjust = 1.0, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 14,
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "grey20", size = 12,
                                    angle = 0, hjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) + labs(tag = "A")+
  guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "top", direction = "verticle"))
  
# back transform x axis


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

st_all_cont_minus9 <- st_all_cont_minus9 %>%
  mutate(log_N = log(Nitrogen+1))
# log linear model 
mSTl_lme_log <- nlme::lme(stability ~  log_N * disk + field,
                          random = (~ 1 | grid), data = st_all_cont_minus9)
# log quadratic model
mSTq_lme_log <- nlme::lme(stability ~  disk * poly(log_N,2,raw=TRUE) + field, 
                          random = (~ 1 | grid), data= st_all_cont_minus9)

# test model fit using AIC function
rawaic <- AIC(mSTl_lme, mSTq_lme, mSTl_lme_log, mSTq_lme_log)
nR <- dim(st_all_cont_minus9)[1]
aictable(rawaic, nR) #  log linear model fit best

MuMIn::r.squaredGLMM(mSTl_lme_log)
an.mSTl_log <- anova(mSTl_lme_log)
summary(mSTl_lme_log)
anova(mSTl_lme_log)

# determine the average trend across fields for plotting purposes
cfa_ST <- ggeffects::ggemmeans(mSTl_lme_log, terms=c("log_N", "disk")) %>%
  as_tibble() %>%
  mutate(xo = expm1(x))


# plot new predicted lines
Fig1B_newmod<- ggplot() +
  geom_point(data = st_all_cont_minus9, aes(x=Nitrogen, y=stability, 
                                            group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = cfa_ST, aes(x = xo, y = predicted, 
                                   group = group, color = group),linewidth = 1) +
  geom_ribbon(data = cfa_ST, aes(
    x = xo,
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
                      labels=c("Intact", "Disturbed"))+
  ylab("Stability")+
  xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
  scale_x_continuous(breaks = c(0, 1, 2, 3.4, 5.4, 9.5, 17.0, 27.2), labels = c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2")) +
  lims(y=c(.5,4.25))+
  labs(legend="Disturbance")+
  theme_bw()+
  #coord_fixed()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, 
                                   angle = 45, hjust = 1.0, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 14, 
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 14, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "grey20", size = 12,
                                    angle = 0, hjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  #theme(legend.position="top") 
  labs(tag = "A")+
  guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "top", direction = "verticle"))



# save ggplot objects to combine plots from another script
saveRDS(Fig1A_newmod, file = here::here("data/syncvsN.rds"))
saveRDS(Fig1B_newmod, file = here::here("data/stabvsN.rds"))




# combine to produce figure as a pdf
# legend_fig1 <- get_legend(Fig1B_newmod)
# 
# quartz(width = 10, height = 5)
# 
# Figure1ab<- Fig1A_newmod / Fig1B_newmod + 
#   plot_layout(ncol=1, nrow= 3, widths = c(8,8), height = c(12,12)) & 
#   theme(legend.position = "none")
# 
# # change plot_layout arguments to adjust heights 
# pdf(file = "Figures/Figure1_AB.pdf",width = 8, height = 12)
# Figure1ab / as_ggplot(legend_fig1) + plot_layout(height=c(12,12,2))
# dev.off()



