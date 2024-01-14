
#############################################################
# FIGURE 2 - SYNCHRONY AND STABILITY BREAK DOWN
#############################################################

source(here::here("data_analyses/Source_figures.R"))

Fig2A<- ggplot()+
  geom_point(data = popcomvar_exp12_2, mapping =aes(x=pop,y=comm, col=Nitrogen,shape=disk),alpha=0.4)+
  geom_point(data= avgpopcommvar, mapping = aes(x=meanpop,y=meancomm, fill=Nitrogen,shape=disk),size=3)+
  #transition_reveal(cvcommip)+
  geom_abline(slope = 1)+
  scale_x_continuous(limits=c(0,0.6))+
  scale_y_continuous(limits=c(0,0.6))+
  labs(x="Aggregate Population \n Variability", y="Community Variability",tag = "A")+
  coord_fixed()+
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = 'bottom')+
  geom_line(data = avgpopcommvar, mapping = aes(x=meanpop, y=meancomm, group=Nitrogen, col=Nitrogen))+
  scale_shape_manual(name = "Disturbance",
                     labels = c("Intact in 1982", "Disturbed in 1982"),
                     values=c(21,24))+
  scale_fill_viridis_d(option = "D",direction=-1, na.value="grey72",
                       name="Nitrogen Addition",
                       breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                       labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))+
  scale_colour_viridis_d(option = "D",direction=-1, na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+ 
  guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "left", direction = "verticle"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"))) 
Fig2A

#### FIGURE 3B - mean over standard devation of total biomass ####
{
  ref_stability <- 161.8789/62.76717
  x <- seq(0,500)
  y <- x*ref_stability
  ref_line <- as.data.frame(cbind(x,y))
  
  ggplot() +
    geom_line(data=ref_line, mapping=aes(x=x, y=y))
  
  Fig2B<- ggplot()+
    geom_point(data = biomass_all_cont_minus9andoutliers, mapping =aes(x=stdev_biomass,y=mean_biomass, col=Nitrogen,shape=disk),alpha=0.4)+
    geom_point(data = biomass_disk_N, mapping =aes(x=meanofsd_biomass,y=meanofmean_biomass, fill=Nitrogen,shape=disk),size=3)+
    #transition_reveal(cvcommip)+
    #geom_abline(slope = 1, intercept = 99.11173)+
    #geom_line(data = ref_line, mapping=aes(x=x, y=y), color="black")+
    geom_abline(slope = ref_stability)+
    #scale_x_continuous(limits=c(50,300))+
    #scale_y_continuous(limits=c(100,400))+
    labs(x="Standard Deviation \n of Total Biomass", y="Mean of Total Biomass",tag = "B")+
    #coord_fixed()+
    theme_bw() +
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
          legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(),
          legend.position = 'bottom')+
    geom_line(data = biomass_disk_N, mapping = aes(x=meanofsd_biomass, y=meanofmean_biomass, group=Nitrogen, col=Nitrogen))+
    scale_shape_manual(name = "Disturbance",
                       labels = c("Intact in 1982", "Disturbed in 1982"),
                       values=c(21,24))+
    scale_fill_viridis_d(option = "D",direction=-1, na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))+
    scale_colour_viridis_d(option = "D",direction=-1, na.value="grey72",
                           name="Nitrogen Addition",
                           breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                           labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))+
    guides(fill=guide_legend(override.aes=list(shape=21)))+ 
    guides(shape = guide_legend(override.aes = list(size = 3), title = "Disturbance", title.position = "left", direction = "verticle"))+ 
    guides(theme(legend.title = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
                 legend.text=element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"))) 
  Fig2B
}
# calculate mean as mean(total_biomass) across time
# calculate stdev as sd(total_biomass) across time
# summarize (mean_of_stab = mean(mean), std_of_stab=mean(std))

Fig2_fin <-  Fig2A + Fig2B + plot_layout(ncol = 2, guides = 'collect') & 
  theme(legend.position = 'bottom') 

pdf(file = "Figures/Figure2_AB.pdf", width = 10, height = 6)
Fig2_fin
dev.off()



# modelling the effect of nitrogen and disk on population variability and community variability
popcomvar_exp12_2$Nitrogen <- as.numeric(popcomvar_exp12_2$Nitrogen)
pop_mod_int <- lm(pop ~ Nitrogen*disk + field, data=popcomvar_exp12_2)
pop_mod_noint <- lm(pop ~ Nitrogen+disk + field, data=popcomvar_exp12_2)
comm_mod_int <- lm(comm ~ Nitrogen*disk + field, data=popcomvar_exp12_2)
comm_mod_noint <- lm(comm ~ Nitrogen+disk + field, data=popcomvar_exp12_2)
summary(pop_mod_int)
summary(pop_mod_noint)
summary(comm_mod_int)
summary(comm_mod_noint)
AIC(pop_mod_int, pop_mod_noint)
AIC(comm_mod_int,comm_mod_noint)
xtable(comm_mod_noint)
xtable(pop_mod_noint)

biomass_all_cont_minus9$Nitrogen <- as.factor(biomass_all_cont_minus9$Nitrogen)
mean_mod <- lm(mean_biomass~Nitrogen*disk + field, data = biomass_all_cont_minus9)
std_mod <- lm(stdev_biomass~Nitrogen*disk + field, data = biomass_all_cont_minus9)
summary(mean_mod)
summary(std_mod)
xtable(std_mod)
stab_mod <- lm(hand_stab~Nitrogen*disk+field, data = biomass_all_cont_minus9)
summary(stab_mod)
xtable(mean_mod)
