
######################################
# FIGURE 5 - biomass through time
######################################

# MD supplemental figure time series 

# New facet label names for disk and ntrt ###
exp.labs <- c("Intact in 1982", "Disturbed in 1982")
names(exp.labs) <- c("1","2")

#ntrt.labs <- c(expression(Control (0 g N/" ~m^2 "/ year)"), expression("High  (9.5 g N/" ~m^2 "/ year"))

ntrt.labs <- c("Control (0 g N)", "High (9.5 g N)")
names(ntrt.labs) <- c("1", "6")

top5sp_avg_subset$year<-as.numeric(top5sp_avg_subset$year)
topbyfunc_subset$year<-as.numeric(topbyfunc_subset$year)
totalbiomass_subset$year<-as.numeric(totalbiomass_subset$year)

biomasstopfunc<-topbyfunc_subset%>% ggplot(aes(x=year, y=meanbiomass, group=species, color=species))+
  geom_line() + 
  geom_line(data=totalbiomass_subset, aes(x = year, y=meantotalbiomass))+
  facet_grid(vars(ntrt),vars(exp), labeller = labeller(ntrt = ntrt.labs, exp = exp.labs))+
  ylab(expression(paste("Mean aboveground biomass (g / ", m^2,")")))+
  xlab("")+
  scale_colour_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4","#00A9FF","#C77CFF","#FF61CC", "black"), labels=c('Agr. rep. (C3 grass)', 'Art. lud. (Pe. forb)','Lath. ven. (Legume)', 'Poa pra. (C3 grass)','Poly. conv. (An. forb)','Rub. sp. (Shrub)', 'Schiz. scop. (C4 grass)', 'Soli. rig. (Pe. forb)', 'Total Biomass'))+
  annotate("rect", xmin = 1982, xmax = 1988, ymin = 0, ymax = 550,alpha = .1,fill = "blue")+
  annotate("rect", xmin = 1994, xmax = 2004, ymin = 0, ymax = 550,alpha = .1,fill = "blue")+
  annotate(geom="text", x=1985, y=500, label="Transient \n phase", size=4)+
  annotate(geom="text", x=1999, y=500, label="Post-transient \n phase", size=4)+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 8,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  theme(legend.position= "bottom")+
  theme(strip.text = element_text(size = 16))+
  theme(legend.margin=margin(0,1,1,0))+
  theme(legend.box = "vertical")+
  theme(legend.box.margin=margin(-10,1,1,-10))+
  guides(shape = guide_legend(override.aes = list(size = 0.5), nrow=3,byrow=TRUE))
biomasstopfunc

# save to repo
pdf(file = "Figures/Figure5_timeseries.pdf", width = 7.5, height = 6)
biomasstopfunc
dev.off()

### log version of above figure##

topbyfunc_subset$meanbiomass1<-1+topbyfunc_subset$meanbiomass

logbiomasstopfunc<-topbyfunc_subset%>% ggplot(aes(x=year, y=meanbiomass1, group=species, color=species))+
  geom_line() + 
  geom_line(data=totalbiomass_subset, aes(x = year, y=meantotalbiomass))+
  facet_grid(vars(ntrt),vars(exp), labeller = labeller(ntrt = ntrt.labs, exp = exp.labs))+
  ylab(expression(paste("Mean aboveground biomass (g / ", m^2,")")))+
  xlab("")+
  scale_y_continuous(trans='log10') +
  scale_colour_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4","#00A9FF","#C77CFF","#FF61CC", "black"), labels=c('Agr. rep. (C3 grass)', 'Art. lud. (Pe. Forb)','Lath. ven. (Legume)', 'Poa pra. (C3 grass)','Poly. conv. (An. forb)','Rub. sp. (Woody)', 'Schiz. scop. (C4 grass)', 'Soli. rig. (Forb)', 'Total Biomass'))+
  annotate("rect", xmin = 1982, xmax = 1988, ymin = 0, ymax = 550,alpha = .1,fill = "blue")+
  annotate("rect", xmin = 1994, xmax = 2004, ymin = 0, ymax = 550,alpha = .1,fill = "blue")+
  annotate(geom="text", x=1985, y=500, label="Transient \n phase", size=4)+
  annotate(geom="text", x=1999, y=500, label="Post-transient \n phase", size=4)+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(legend.position = "right")+theme(legend.title= element_blank())+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_blank(),
        legend.text = element_text(color = "grey20", size = 8,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
  theme(legend.position= "bottom")+
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(size = 16))+
  theme(legend.margin=margin(0,1,1,0))+
  theme(legend.box = "vertical")+
  theme(legend.box.margin=margin(-10,1,1,-10))+
  guides(shape = guide_legend(override.aes = list(size = 0.5), nrow=3,byrow=TRUE))

# save to repo
pdf(file = "Figures/Figure5_timeseriesLOGSCALE.pdf", width = 7.5, height = 6)
logbiomasstopfunc
dev.off()