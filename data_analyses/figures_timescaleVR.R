# Figures ####
# Figures from Zhao et al 2020 - Both experiments - All nutrients 

source(here('VR_calculations.R'))

# paired plots from short to long timescales, each response by nutrient group & experiment:
timescale_comm <- ggplot(TSVR_all_wide)+
  geom_point(aes(x=timescale, y=as.numeric(CV_comm),color = timescale), size = 2)+
  geom_line(aes(x=timescale, y=as.numeric(CV_comm), group = uniqueID))+
  theme(legend.position = "none")+
  ylab("community variability (CV2com)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) +
  facet_grid(cols = vars(ntrt_group), rows = vars(exp)) + 
  theme_bw() +
  theme(legend.position = 'none')+
  scale_color_manual(values = c('turquoise3','tomato1'))

timescale_ip <- ggplot(TSVR_all_wide)+
  geom_point(aes(x=timescale, y=as.numeric(CV_ip),color = timescale), size = 2)+
  geom_line(aes(x=timescale, y=as.numeric(CV_ip), group = uniqueID))+
  theme(legend.position = "none")+
  ylab("aggregate population variability (CV2com_ip)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  facet_grid(cols = vars(ntrt_group), rows = vars(exp)) + 
  theme_bw() +
  theme(legend.position = 'none')+
  scale_color_manual(values = c('turquoise3','tomato1'))

weightedvr <- ggplot(TSVR_all_wide)+
  geom_point(aes(x=timescale, y=as.numeric(weighted),color = timescale), size = 2)+
  geom_line(aes(x=timescale, y=as.numeric(weighted), group = uniqueID))+
  theme(legend.position = "none")+
  ylab("timescale-specific variance ratio")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  facet_grid(cols = vars(ntrt_group), rows = vars(exp)) + 
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('turquoise3','tomato1'))

# figure CV pop vs CV comm with nutrient groups & disturbance
CV_CV_short <- ggplot(filter(TSVR_all_wide, timescale %in%'short'))+
  geom_point(aes(x=as.numeric(CV_ip), y=as.numeric(CV_comm)), size = 1.5, col = 'turquoise3')+
  geom_abline(intercept = 0, slope = 1)+
  xlab("")+
  ylab("community variability (CV2com)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  facet_grid(cols = vars(ntrt_group), rows = vars(exp)) + 
  theme_bw()
CV_CV_long <-ggplot(filter(TSVR_all_wide, timescale %in%'long'))+
  geom_point(aes(x=as.numeric(CV_ip), y=as.numeric(CV_comm)), size = 1.5, col = 'tomato1')+
  geom_abline(intercept = 0, slope = 1)+
  xlab("aggregate population variability (CV2com_ip)")+
  ylab("community variability (CV2com)")+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  facet_grid(cols = vars(ntrt_group), rows = vars(exp)) + 
  theme_bw()
library(patchwork)
timescale_VR_fig <- CV_CV_short / CV_CV_long
