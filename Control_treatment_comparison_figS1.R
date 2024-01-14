#####################################
# This script produces figure S1 to compare two control treatments, one with 
# micronutrients added, and one without micronutrients added
########################

# Read in data and functions from source code
source(here::here("data_cleaning/subsetting_CC.R"))

#### Fig S1A. Comparison of control conditions on synchrony ####
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


# compare synchrony for the two control treatments 
VR_all_cont_control <- VR_all_cont %>% 
  dplyr::filter(ntrt %in% c(1,9))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

Fig_control_VR <- ggplot(data=VR_all_cont_control, 
                         aes(x=ntrt, y=VR, fill = disk)) +
  geom_violin(position=position_dodge(1)) +
  stat_summary(fun.data=data_summary, position=position_dodge(1)) +
  scale_fill_manual(values = c("#D55E00", "skyblue"),
                    name="Disturbance",
                    breaks=c("0", "1"),
                    labels=c("Intact in 1982", "Disturbed in 1982"))+
  ylab("Synchrony")+
  xlab("") + 
  geom_hline(yintercept=1, color="darkgrey", linetype = "dashed") +
  theme_bw()+
  lims(y=c(0,1.6))+
  theme(axis.text.x = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, vjust = 0, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 12, 
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5,2.0)) +   
  scale_x_discrete(labels=c(expression("0 N + "~mu, "0 N + 0"~mu))) +
  labs(tag = "A")

#### Fig S1B. Comparison of control conditions on stability ####
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

st_all_cont_control <- st_all_cont %>% 
  dplyr::filter(ntrt %in% c(1,9))

Fig_control_stability <- ggplot(data=st_all_cont_control, 
                                aes(x=ntrt, y=stability, fill = disk)) +
  geom_violin(position=position_dodge(1)) +
  stat_summary(fun.data=data_summary, position=position_dodge(1)) +
  scale_fill_manual(values = c("#D55E00", "skyblue"),
                    name="Disturbance",
                    breaks=c("0", "1"),
                    labels=c("Intact in 1982", "Disturbed in 1982"))+
  ylab("Stability")+
  xlab("Treatment") +
  theme_bw()+
  lims(y=c(.5,4.25)) +
  theme(axis.text.x = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   angle = 0, hjust = .5, vjust = 0, 
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 12, 
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) +
  scale_x_discrete(labels=c(expression("0 N + "~mu, "0 N + 0"~mu))) + 
  labs(tag = "B")

#Create final figure
legend_controlfig <- get_legend(Fig_control_stability)
quartz(width = 6, height = 4)
control_plots<- Fig_control_VR / Fig_control_stability + 
  plot_layout(ncol=1, widths = c(2,2)) & theme(legend.position = "none")

#Save figure as pdf
pdf(file = "Figures/Controlfigures_AB.pdf",width = 4, height = 6)
control_plots / as_ggplot(legend_controlfig) + plot_layout(heights=c(2,2,2))
dev.off()

#Linear model comparing micronutrient effect on synchrony
VR_control_model <- lm(VR ~  ntrt*disk + field, data=VR_all_cont_control)
summary(VR_control_model)
plot(VR_control_model)

#Linear model comparing micronutrient effect on stability
stability_control_model <- lm(stability ~  ntrt*disk + field, 
                              data=st_all_cont_control)
summary(stability_control_model)
plot(stability_control_model)