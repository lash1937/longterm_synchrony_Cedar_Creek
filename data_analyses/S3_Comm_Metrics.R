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
