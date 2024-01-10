#####################################
# This script produces all supplemental figures.
#####################################

#Read in data, functions and packages
source(here::here("data_analyses/Source_figures.R"))


#### Fig S1. Comparison of control conditions ####
#Figure comparing micronutrient effects on synchrony
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

#Figure comparing micronutrient effects on stability
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

#### Fig S2. #### 

#### Fig S3. Comparison of community metrics across succession ####
#prep dataframe for figure S3 creation
sem.a.df <- SEM.a.df  
sem.a.df$transience <- c("transient")
sem.b.df <- SEM.b.df
sem.b.df$transience <- c("post-transient")
all_bothtime <- rbind(sem.a.df, sem.b.df)
all_bothtime$transience <- as.factor(all_bothtime$transience)
levels(all_bothtime$transience) <- c("Transient", "Post-Transient")

## Richness
#Filter dataset for richness
all_bothtime_r <- all_bothtime %>%
  select(Nitrogen, Richness,Disturbance, transience)
#Plot effects of global change on richness across succession
FigS3Richness_bothtime<- ggplot() +
  geom_point(data =all_bothtime_r, 
             aes(x=Nitrogen, y=Richness, group = Disturbance, 
                 col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_r, 
              aes(x=Nitrogen, y=Richness, group = Disturbance, 
                  col = Disturbance),method="lm",se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 10, 
                                   angle = 0, hjust = 0, face = "plain"),
        legend.position = c(0.77, 0.80),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  labs(x = "", y = "Species Richness", tag = "A") +
  theme(strip.text.x = element_text(size = 16))

## Evenness
#Filter dataset for evenness
all_bothtime_e <- all_bothtime %>%
  select(Nitrogen, Evenness,Disturbance, transience)
#Plot effects of global change on evenness across succession
FigS3Evenness_bothtime<- ggplot() +
  geom_point(data =all_bothtime_e, 
             aes(x=Nitrogen, y=Evenness, group = Disturbance, 
                 col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_e, 
              aes(x=Nitrogen, y=Evenness, group = Disturbance, 
                  col = Disturbance), method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982")) +
  facet_wrap(~transience) +
  theme_bw() +
  theme(axis.title.y = element_text(color = "grey20", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none", 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  labs(x = "", y = "Species Evenness", tag = "B") +
  theme(strip.text.x = element_blank())

## Synchrony
#Filter dataset for VR
all_bothtime_vr <- all_bothtime %>%
  select(Nitrogen,VR,Disturbance, transience)
#Plot effects of global change on synchrony across succession
FigS3Synchrony_bothtime<- ggplot() +
  geom_point(data =all_bothtime_vr, 
             aes(x=Nitrogen, y=VR, group = Disturbance, 
                 col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_vr, 
              aes(x=Nitrogen, y=VR, group = Disturbance, 
                  col = Disturbance), method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none", 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  labs(x = "", y = "Synchrony", tag = "C") +
  theme(strip.text.x = element_blank())

## Stability
#Filter dataset for stability
all_bothtime_stab <- all_bothtime %>%
  select(Nitrogen,Stability,Disturbance, transience)
#Plot effects of global change on synchrony across succession
FigS3Stability_bothtime<- ggplot() +
  geom_point(data =all_bothtime_stab, 
             aes(x=Nitrogen, y=Stability, group = Disturbance, 
                 col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_stab, 
              aes(x=Nitrogen, y=Stability, group = Disturbance, 
                  col = Disturbance),method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, 
                                    angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none", 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_text(color = "grey20", size = 12, 
                                  angle = 0, hjust = .5, face = "plain"),
        axis.text.x=element_text(color = "grey20", size = 10, 
                                 angle = 0, hjust = .5, face = "plain")) +
  labs(x = "", y = "", tag = "D") +
  ylab("Stability")+
  xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_blank())

#create pdf of final figure S3
pdf(file = "Figures/SupFigureS3_ABCD_bothtime.pdf",width = 10, height = 10)
FigS3Richness_bothtime / FigS3Evenness_bothtime / FigS3Synchrony_bothtime / FigS3Stability_bothtime
dev.off()
