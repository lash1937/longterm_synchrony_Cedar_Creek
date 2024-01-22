
#####################################
# This script produces figure 2, breaking down total biomass across functional groups
#####################################

# Read in data and functions from source code
source(here::here("data_cleaning/subsetting_CC.R"))

# prepare dataframe for figure
unique_ID_long2 <- unique_ID_long %>%  tidyr::separate("uniqueID", c("field", "exp", "plot", "disk", "ntrt"), sep="_") %>%  separate("expntrtyear", c("exp2", "ntrt2", "year"), sep="_") %>% select(field, exp, ntrt, plot, year, Species, Abundance )

# remove extra text
unique_ID_long2$species<-str_remove(unique_ID_long2$Species, "mass.above.")

# summarize biomass all to get most abundance spp 
biomassbysp<-unique_ID_long2 %>% dplyr::group_by(species) %>% dplyr::summarise(totalbiomass=sum(Abundance, na.rm=TRUE))
biomassbysp2 <- merge(biomassbysp,sp.df, by='species', all.x=TRUE)

# pull out the top 5 most abundance species
top5sp<-subset(unique_ID_long2, species=="AGROPYRON REPENS"|Species=="POA PRATENSIS"|Species=="ARTEMISIA LUDOVICIANA"|Species=="SCHIZACHYRIUM SCOPARIUM"|Species=="SOLIDAGO RIGIDA")
topbyfunc<-subset(unique_ID_long2, species=="AGROPYRON REPENS"|species=="POA PRATENSIS"|species=="ARTEMISIA LUDOVICIANA"|species=="SCHIZACHYRIUM SCOPARIUM"|species=="SOLIDAGO RIGIDA"|species=="RUBUS SP."|species=="LATHYRUS VENOSUS"|species=="POLYGONUM CONVOLVULUS")
totalbiomass<-unique_ID_long2 %>% dplyr::group_by(field, exp,year, plot, ntrt) %>% dplyr::summarise(totalbiomass=sum(Abundance, na.rm=TRUE))

# average across treatments
top5sp_avg<-top5sp%>% dplyr::group_by(exp,year,ntrt, species) %>% dplyr::summarise(meanbiomass=mean(Abundance))
topbyfunc_avg<-topbyfunc%>% dplyr::group_by(exp,year,ntrt, species) %>% dplyr::summarise(meanbiomass=mean(Abundance))
totalbiomass_avg<-totalbiomass %>% dplyr::group_by(exp,year,ntrt) %>% dplyr::summarise(meantotalbiomass=mean(totalbiomass))

# subset to only control and high nutrient 
top5sp_avg_subset<-subset(top5sp_avg, ntrt==1|ntrt==6)
topbyfunc_subset<-subset<-subset(topbyfunc_avg, ntrt==1|ntrt==6)
top5sp_avg_subset$ntrt <- factor(top5sp_avg_subset$ntrt, levels = c("1", "6"))
topbyfunc_subset$ntrt <- factor(topbyfunc_subset$ntrt, levels = c("1", "6"))
totalbiomass_subset<-subset(totalbiomass_avg, ntrt==1|ntrt==6)
totalbiomass_subset$ntrt <- factor(totalbiomass_subset$ntrt, levels = c("1", "6"))
totalbiomass_subset<-totalbiomass_subset %>% dplyr::mutate(species="total biomass")

# New facet label names for disk and ntrt ###
exp.labs <- c("Intact in 1982", "Disturbed in 1982")
names(exp.labs) <- c("1","2")
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
pdf(file = "Figures/Figure2_timeseries.pdf", width = 7.5, height = 6)
biomasstopfunc
dev.off()

