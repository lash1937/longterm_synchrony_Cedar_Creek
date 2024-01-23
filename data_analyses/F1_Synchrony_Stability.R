
#####################################
# This script produces summary outputs of linear and quadratic models,
# which are tested to fit the relationships between global change drivers
# and synchrony and stability. The figures produced represent the effects 
# of nitrogen addition and soil disturbance on synchrony and stability.
#####################################

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
VR_all_cont$Nitrogen<-as.numeric(as.character(VR_all_cont$Nitrogen))

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
st_all_cont$Nitrogen<-as.numeric(as.character(st_all_cont$Nitrogen))

# subset out nutrient treatment 9, control plots without micronutrients
st_all_cont_minus9 <- st_all_cont%>%
  dplyr::filter(ntrt < 9)

# unite unique ID identifers in new dataset
st_all_cont_minus9 <- st_all_cont_minus9 %>% 
  tidyr::unite("uniqueID", 1:5, sep="_",remove = FALSE)

# soil disturbance treatment, disk, as a factor
st_all_cont_minus9$disk <- as.factor(st_all_cont_minus9$disk)

#### ISSUE WITH MODEL PREDICTION ----------------------------------------------
#### Fig 1A. The effect of global change drivers on synchrony ####
# compare a linear model to a quadratic model to best describe the relationship 
# between synchrony and global change drivers

# create unique grid variable
VR_all_cont_minus9 <- VR_all_cont_minus9 %>%
  mutate(
    grid = factor(paste0(field, exp)))

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

# refit the model with contrast / sum / deviations coding for field
# so that the plot will show the curve for the "average field"
VR_all_cont_minus9$field_contr <- as.factor(VR_all_cont_minus9$field)
contrasts(VR_all_cont_minus9$field_contr) <- contr.sum(
  length(levels(VR_all_cont_minus9$field_contr)))

mVRl_contr <- nlme::lme(VR ~  Nitrogen * disk + field,
                        random = (~ 1 | grid), data = VR_all_cont_minus9)
        
# create new model matrix for plotting fitted line
n_per_disk <- 100
Xmm_pred <- data_frame(
Intercept = rep(1, n_per_disk * 2),
disk1 = rep(c(0, 1), each = n_per_disk),
N1 = rep(seq(0, 30, length.out = n_per_disk), 2),
field1 = rep(0, n_per_disk * 2),
field2 = field1,
disk_N1 = disk1 * N1)

Xmm_pred <- as.matrix(Xmm_pred)
V_mVRl_contr <- vcov(mVRl_contr)
 
beta_hat_vr <- mVRl_contr$coefficients
beta_hat_vr <- beta_hat_vr$fixed
t_star_vr <- qt(0.975, df = 208)

df_pred_vr <- tibble(
y_pred = as.double(Xmm_pred %*% beta_hat_vr),
se_y = as.double(sqrt(diag(Xmm_pred %*% V_mVRl_contr %*% t(Xmm_pred)))),
low = y_pred - t_star_vr * se_y,
high = y_pred + t_star_vr * se_y,
Nitrogen = Xmm_pred[, "N1"],
disk = as.factor(Xmm_pred[, "disk1"]))

# find the minimum
VR_predict_min <- df_pred_vr %>%
 dplyr::group_by(disk) %>%
 slice(which.min(y_pred))

# # plot new predicted lines to smooth the quadratic
Fig1A_newmod<- ggplot() +
  geom_point(data = VR_all_cont_minus9, aes(x=Nitrogen, y=VR, group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = df_pred, aes(x = Nitrogen, y = y_pred, group = disk, 
                                color = disk),linewidth = 1) +
  geom_ribbon(data = df_pred, aes(
    x = Nitrogen,
    y = y_pred,
    group= disk,
    fill = disk,
    ymin = low,
    ymax = high),
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
  mutate(
    grid = factor(paste0(field, exp))
  )

# linear model
mSTl_lme <-
  nlme::lme(stability ~ disk * Nitrogen + field, random = (~1 | grid), data= st_all_cont_minus9)
# quadratic model
mSTq_lme <-  nlme::lme(stability ~ disk * poly(Nitrogen,2,raw=TRUE) + field,random = (~1 | grid), data = st_all_cont_minus9)

# test model fit using AIC function
rawaic <- AIC(mSTl_lme, mSTq_lme)
nR <- dim(st_all_cont_minus9)[1]
aictable(rawaic, nR) # linear model fit best

# refit the model with contrast / sum / deviations coding for field
#  so that the plot will show the curve for the "average field"
st_all_cont_minus9$field_contr <- as.factor(st_all_cont_minus9$field)
contrasts(st_all_cont_minus9$field_contr) <- contr.sum(
  length(levels(st_all_cont_minus9$field_contr)))

mSTl_contr <- nlme::lme(stability ~ disk * Nitrogen + field, random = (~1 | grid), data= st_all_cont_minus9)

n_per_disk <- 100
Xmm_pred <- data_frame(
  Intercept = rep(1, n_per_disk * 2),
  disk1 = rep(c(0, 1), each = n_per_disk),
  N1 = rep(seq(0, 30, length.out = n_per_disk), 2),
  field1 = rep(0, n_per_disk * 2),
  field2 = field1,
  disk_N1 = disk1 * N1)

Xmm_pred <- as.matrix(Xmm_pred)
V_mSTl_contr <- vcov(mSTl_contr)
beta_hat_st <- mSTl_contr$coefficients
beta_hat_st <- beta_hat_st$fixed
t_star_st <- qt(0.975, df = 210)

df_pred_st <- tibble(
  y_pred = as.double(Xmm_pred %*% beta_hat_st),
  se_y = as.double(
    sqrt(diag(Xmm_pred %*% V_mSTl_contr %*% t(Xmm_pred)))),
  low = y_pred - t_star_st * se_y,
  high = y_pred + t_star_st * se_y,
  Nitrogen = Xmm_pred[, "N1"],
  disk = as.factor(Xmm_pred[, "disk1"]))

# find the minimum
st_predict_min <- df_pred_st %>%
  dplyr::group_by(disk) %>%
  slice(which.min(y_pred))

# plot new predicted lines
Fig1B_newmod<- ggplot() +
  geom_point(data = st_all_cont_minus9, aes(x=Nitrogen, y=stability, 
                                            group = disk, 
                                            col = disk), shape = 21) +
  geom_line(data = df_pred_st, aes(x = Nitrogen, y = y_pred, 
                                group = disk, color = disk),linewidth = 1) +
  geom_ribbon(data = df_pred_st, aes(
    x = Nitrogen,
    y = y_pred,
    group= disk,
    fill = disk,
    ymin = low,
    ymax = high),
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
quartz(width = 10, height = 5)
library(patchwork)

Figure1ab<- Fig1A_newmod / Fig1B_newmod + 
  plot_layout(ncol=1, widths = c(3,3)) & 
  theme(legend.position = "none")

pdf(file = "Figures/Figure1_AB.pdf",width = 4, height = 6)
Figure1ab / as_ggplot(legend_fig1) + plot_layout(heights=c(2,2,2))
dev.off()


# produce final summary tables for both models 
xtable(summary(mVRq_contr))
xtable(summary(mSTl_contr))
