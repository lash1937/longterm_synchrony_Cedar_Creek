
#######################
# SUPPLEMENTAL FIGURES
#######################

source(here::here("data_analyses/Source_figures.R"))

###########################
### rank abundance curve##
###########################

biomassbysp3<-biomassbysp2  %>% arrange(desc(totalbiomass)) %>% dplyr::filter(totalbiomass > 0)%>% mutate(totalbiomass2=totalbiomass+1) %>% mutate(spyesno = ifelse(str_detect(species, fixed("SP.")), "yes", "no"))

write.csv(biomassbysp3, here("data/spabundance.csv"), row.names=FALSE)

options(scipen = 999)

rankgraph<-ggplot(data=biomassbysp3, aes(x = reorder(species, -totalbiomass), y = totalbiomass2, fill=spyesno)) +
  geom_bar(stat="identity",width = 0.85) +
  scale_y_continuous(trans='log10')+xlab("")+ylab("Log Total biomass") + 
  theme_bw()+  scale_fill_manual(values=c("#00BFC4", "#F8766D"),name="Genus level") +
  theme(axis.text.x = element_text(color = "grey20", size = 7, angle = 90, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        axis.ticks.x = element_blank(),
        #legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 10,angle = 0, hjust = 0, face = "plain"),
        legend.position = c(0.8, 0.8),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())

# save to repo
pdf(file = "Figures/rankabundance.pdf", width = 15, height = 6)
rankgraph
dev.off()

###################################
# ALL DYNAMICS PANELS 2-time
###################################

# KM recreate figure S3
sem.a.df <- SEM.a.df  
sem.a.df$transience <- c("transient")
sem.b.df <- SEM.b.df
sem.b.df$transience <- c("post-transient")
all_bothtime <- rbind(sem.a.df, sem.b.df)

all_bothtime$transience <- as.factor(all_bothtime$transience)
levels(all_bothtime$transience) <- c("Transient", "Post-Transient")

# Richness
all_bothtime_r <- all_bothtime %>%
  select(Nitrogen, Richness,Disturbance, transience)

FigS3Richness_bothtime<- ggplot() +
  geom_point(data =all_bothtime_r, aes(x=Nitrogen, y=Richness, group = Disturbance, col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_r, aes(x=Nitrogen, y=Richness, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  #xlab("Nitrogen")+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 10,angle = 0, hjust = 0, face = "plain"),
        legend.position = c(0.77, 0.80),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  labs(x = "", y = "Species Richness", tag = "A") +
  theme(strip.text.x = element_text(size = 16))

# Evenness
all_bothtime_e <- all_bothtime %>%
  select(Nitrogen, Evenness,Disturbance, transience)

FigS3Evenness_bothtime<- ggplot() +
  geom_point(data =all_bothtime_e, aes(x=Nitrogen, y=Evenness, group = Disturbance, col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_e, aes(x=Nitrogen, y=Evenness, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  #xlab("Nitrogen")+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, face = "plain"),
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

# Synchrony
all_bothtime_vr <- all_bothtime %>%
  select(Nitrogen,VR,Disturbance, transience)

FigS3Synchrony_bothtime<- ggplot() +
  geom_point(data =all_bothtime_vr, aes(x=Nitrogen, y=VR, group = Disturbance, col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_vr, aes(x=Nitrogen, y=VR, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  #xlab("Nitrogen")+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, face = "plain"),
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

# Stability
all_bothtime_stab <- all_bothtime %>%
  select(Nitrogen,Stability,Disturbance, transience)

FigS3Stability_bothtime<- ggplot() +
  geom_point(data =all_bothtime_stab, aes(x=Nitrogen, y=Stability, group = Disturbance, col = Disturbance), shape = 21) +
  geom_smooth(data =all_bothtime_stab, aes(x=Nitrogen, y=Stability, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
  scale_colour_manual(values = c("#D55E00", "skyblue"),
                      name="Disturbance",
                      breaks=c("0", "1"),
                      labels=c("Intact in 1982", "Disturbed in 1982"))+
  #xlab("Nitrogen")+
  facet_wrap(~transience)+
  theme_bw()+
  theme(axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none", 
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.x=element_text(color = "grey20", size = 12, angle = 0, hjust = .5, face = "plain"),
        axis.text.x=element_text(color = "grey20", size = 10, angle = 0, hjust = .5, face = "plain")) +
  labs(x = "", y = "", tag = "D") +
  ylab("Stability")+
  xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_blank())

library(patchwork)
pdf(file = "Figures/SupFigureS3_ABCD_bothtime.pdf",width = 10, height = 10)
FigS3Richness_bothtime / FigS3Evenness_bothtime / FigS3Synchrony_bothtime / FigS3Stability_bothtime
dev.off()

#MD original figure S3
### Richness panel 1 

# Fig1Richness_bothtime<- ggplot() +
#   geom_point(data =all_bothtime, aes(x=Nitrogen, y=Richness, group = Disturbance, col = Disturbance), shape = 21) +
#   geom_smooth(data =all_bothtime, aes(x=Nitrogen, y=Richness, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
#   scale_colour_manual(values = c("#D55E00", "skyblue"),
#                       name="Disturbance",
#                       breaks=c("0", "1"),
#                       labels=c("Intact in 1982", "Disturbed in 1982"))+
#   #xlab("Nitrogen")+
#   facet_wrap(~transience)+
#   theme_bw()+
#   theme(axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, face = "plain"),
#         legend.title = element_blank(),
#         legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
#         legend.position = c(0.88, 0.80),
#         panel.grid.minor.y=element_blank(),
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank()) +
#   labs(x = "", y = "Species Richness", tag = "A") +
#   theme(strip.text.x = element_text(size = 16))

####################
### Evenness panel 2
# Fig1Evenness_bothtime<- ggplot() +
#   geom_point(data =all_bothtime, aes(x=Nitrogen, y=Evenness, group = Disturbance, col = Disturbance), shape = 21) +
#   geom_smooth(data =all_bothtime, aes(x=Nitrogen, y=Evenness, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
#   scale_colour_manual(values = c("#D55E00", "skyblue"),
#                       name="Disturbance",
#                       breaks=c("0", "1"),
#                       labels=c("Intact in 1982", "Disturbed in 1982"))+
#   #xlab("Nitrogen")+
#   facet_wrap(~transience)+
#   theme_bw()+
#   theme(axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, face = "plain"),
#         legend.title = element_blank(),
#         legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
#         legend.position = 'none',
#         panel.grid.minor.y=element_blank(),
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank()) +
#   labs(x = "", y = "Species Evenness", tag = "B") +
#   theme(strip.background = element_blank()) +
#   theme(strip.text.x = element_blank())
# Fig1Evenness_bothtime
# 
# ####################
# ### VR panel 3 
# Fig1VR_bothtime<- ggplot() +
#   geom_point(data =all_bothtime, aes(x=Nitrogen, y=VR, group = Disturbance, col = Disturbance), shape = 21) +
#   geom_smooth(data =all_bothtime, aes(x=Nitrogen, y=VR, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
#   scale_colour_manual(values = c("#D55E00", "skyblue"),
#                       name="Disturbance",
#                       breaks=c("0", "1"),
#                       labels=c("Intact in 1982", "Disturbed in 1982"))+
#   scale_y_continuous(breaks = c(0.0,0.5,1.0,1.5,2.0))+
#   ylab("Synchrony")+
#   #xlab("Nitrogen")+
#   geom_hline(yintercept=1, color="darkgrey", linetype = "dashed") +
#   facet_wrap(~transience)+
#   theme_bw()+
#   theme(axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, face = "plain"),
#         legend.title = element_blank(),
#         legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
#         legend.position = 'none',
#         panel.grid.minor.y=element_blank(),
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank()) +
#   labs(x = "", y = "Synchrony", tag = "C") +
#   theme(strip.background = element_blank()) +
#   theme(strip.text.x = element_blank())
# 
# ####################
# ### Stability panel 4 
# Fig1Stability_bothtime<- ggplot() +
#   geom_point(data =all_bothtime, aes(x=Nitrogen, y=Stability, group = Disturbance, col = Disturbance), shape = 21) +
#   geom_smooth(data =all_bothtime, aes(x=Nitrogen, y=Stability, group = Disturbance, col = Disturbance),method="lm", se=FALSE)+
#   scale_colour_manual(values = c("#D55E00", "skyblue"),
#                       name="Disturbance",
#                       breaks=c("0", "1"),
#                       labels=c("Intact in 1982", "Disturbed in 1982"))+
#   facet_wrap(~transience)+
#   theme_bw()+
#   theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, face = "plain"),
#         axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#         axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
#         axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, face = "plain"),
#         legend.title = element_blank(),
#         legend.position = 'none', 
#         panel.grid.minor.y=element_blank(),
#         panel.grid.major.y=element_blank(),
#         panel.grid.minor.x=element_blank(),
#         panel.grid.major.x=element_blank()) +
#   labs(x = "", y = "", tag = "D") +
#   ylab("Stability")+
#   xlab(expression(paste("Nitrogen (g/",m^2,"/year)")))+
#   theme(strip.background = element_blank()) +
#   theme(strip.text.x = element_blank())
# 
# library(patchwork)
# pdf(file = "Figures/SupFigure1_ABCD_bothtime.pdf",width = 10, height = 10)
# Fig1Richness_bothtime / Fig1Evenness_bothtime / Fig1VR_bothtime / Fig1Stability_bothtime
# dev.off()

#################################
# SUPPLEMENTAL
# Slopes and intercepts 
# for stability v synchrony
#################################


# See Clean_figure 3.R Script