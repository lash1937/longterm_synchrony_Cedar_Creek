#####################################
# This script produces figure 4 that compares the synchrony stability relationship
# across successional timescales and supplemental figure 2 that breaks down the
# slopes and intercepts of this relationship
########################

# Read in data and functions from source code
source(here::here("data_cleaning/subsetting_CC.R"))

#### Fig 4. synchrony stability relationships ####
# prep data for analyses
SEM.b.df$transience<-"Transient"
SEM.a.df$transience<-"Post-transient"
all_bothtime<-rbind(SEM.b.df, SEM.a.df)
all_bothtime$transience<-factor(all_bothtime$transience, levels=c("Transient", "Post-transient"))
all_bothtime$Disturbance<-factor(all_bothtime$Disturbance)

vr_st_df_bothtime <- all_bothtime

vr_st_df_bothtime$Nitrogen <- as.factor(vr_st_df_bothtime$Nitrogen)
vr_st_df_bothtime$Disturbance <- as.factor(vr_st_df_bothtime$Disturbance)
vr_st_df_bothtime$field <- as.factor(vr_st_df_bothtime$field)

# Slopes and intercepts for stability v synchrony
# Loop to find linear models fits of each treatment, for synchrony against stability
vr_st_df_bothtime <- as.data.frame(vr_st_df_bothtime)
library(lme4)
library(ggeffects)

dist <- c(levels(vr_st_df_bothtime$Disturbance))
N <- c(levels(vr_st_df_bothtime$Nitrogen))
trans <- c(levels(vr_st_df_bothtime$transience))
lm_results_bothtime <- matrix(ncol = 9, nrow = 32)
colnames(lm_results_bothtime) <-  c("Disturbance", "transience", "Nitrogen","Slope","CI_lwr_slope","CI_upper_slope", "Intercept", "CI_lwr_int", "CI_upper_int")
k <- 1

for(d in dist) {
  for(t in trans) {
    for(n in N){
      temp <- dplyr::filter(vr_st_df_bothtime, 
                            Disturbance %in% d &
                              transience %in% t &
                              Nitrogen %in% n)
      mod = lm(Stability ~ VR + field, data = temp) 
      em.mod = emtrends(mod, "field", var = "VR")
      em.field = ggemmeans(mod, terms = "VR[0]") # find field intercepts
      lm_results_bothtime[k,] <- c(
        d,t,n,
        as.data.frame(em.mod)[1,2], # Slope
        as.data.frame(em.mod)[1,5], # CI lower slope
        as.data.frame(em.mod)[1,6], # Ci upper slope
        as.data.frame(em.field)[1,2], # intercept
        as.data.frame(em.field)[1,4], # CI lower intercept
        as.data.frame(em.field)[1,5]) # CI upper intercept
      
      k = k + 1
    }
  }
}


lm_results_bothtime <- as.data.frame(lm_results_bothtime)
lm_results_bothtime$Disturbance <- as.factor(lm_results_bothtime$Disturbance)
lm_results_bothtime$transience <- as.factor(lm_results_bothtime$transience)
lm_results_bothtime$Nitrogen <- as.factor(lm_results_bothtime$Nitrogen)
lm_results_bothtime$Slope <- as.numeric(lm_results_bothtime$Slope)
lm_results_bothtime$CI_lwr_slope <- as.numeric(lm_results_bothtime$CI_lwr_slope)
lm_results_bothtime$CI_upper_slope <- as.numeric(lm_results_bothtime$CI_upper_slope)
lm_results_bothtime$Intercept <- as.numeric(lm_results_bothtime$Intercept)
lm_results_bothtime$CI_lwr_int <- as.numeric(lm_results_bothtime$CI_lwr_int)
lm_results_bothtime$CI_upper_int <- as.numeric(lm_results_bothtime$CI_upper_int)

# nitrogen order
levels(lm_results_bothtime$Nitrogen)
lm_results_bothtime$Nitrogen <- factor(lm_results_bothtime$Nitrogen, levels = 
                                         c("0","1","2","3.4","5.4","9.5","17","27.2"))

# transient order
lm_results_bothtime$transience <- factor(lm_results_bothtime$transience, 
                                         levels = c("Transient", "Post-transient"))

# find the verticle lines 'baselines'
baseline <- lm_results_bothtime %>%
  group_by(Disturbance, transience) %>%
  filter(Nitrogen == '0')

# Call Slopes and interceps when N addition = 0
base_slope_intact_trans <- as.numeric(baseline[1,4])
base_int_intact_trans <- as.numeric(baseline[1,7])

base_slope_intact_post <- as.numeric(baseline[2,4])
base_int_intact_post <- as.numeric(baseline[2,7])

base_slope_dist_trans <- as.numeric(baseline[3,4])
base_int_dist_trans <- as.numeric(baseline[3,7])

base_slope_dist_post <- as.numeric(baseline[4,4])
base_int_dist_post <- as.numeric(baseline[4,7])

# rearrange so I can put in dataframe for plotting
baseline_slopes <- c(rep(base_slope_intact_trans, times = 8),
                     rep(base_slope_intact_post, times = 8),
                     rep(base_slope_dist_trans, times = 8),
                     rep(base_slope_dist_post, times = 8))

baseline_ints <- c(rep(base_int_intact_trans, times = 8),
                   rep(base_int_intact_post, times = 8),
                   rep(base_int_dist_trans, times = 8),
                   rep(base_int_dist_post, times = 8))

# attach to dataframe
lm_results_bothtime$baseline_slopes <- NULL
lm_results_bothtime$baseline_ints <- NULL

lm_results_bothtime$baseline_slopes <- baseline_slopes
lm_results_bothtime$baseline_ints <- baseline_ints
tail(lm_results_bothtime)

# Now use slopes and intercepts to create predicted values across VR (spans 0 - 3)
VR_seq <- seq(0,3, length.out = 20)

# create another loop that applies the VR_seq as x for each on the of the model_lines rows
model_response <- matrix(ncol = 12, nrow = 32*20) # 32 unique categoires, for 20 rows of my VR_seq
colnames(model_response) <- c("Disturbance", "transience", "Nitrogen", "Slope", "CI_lwr_slope", "CI_upper_slope", "Intercept",  "CI_lwr_int", "CI_upper_int", "baseline_slopes", "baseline_ints", "predicted") # make a new dataframe where rows are 20 entires for each VR
lm_results_bothtime_mat <- as.matrix(lm_results_bothtime)
dimnames(lm_results_bothtime_mat) <- NULL
k = 1

for(i in 1:32){ # 32 lines
  print(i)
  print(k)
  model_response[as.numeric(k):as.numeric(k+19),1] <- rep(lm_results_bothtime_mat[i,1], times = 20) # pull out the categories for model_lines row (Disturbance, N, Transience)
  model_response[as.numeric(k):as.numeric(k+19),2] <- rep(lm_results_bothtime_mat[i,2], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),3] <- rep(lm_results_bothtime_mat[i,3], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),4] <- rep(lm_results_bothtime_mat[i,4], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),5] <- rep(lm_results_bothtime_mat[i,5], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),6] <- rep(lm_results_bothtime_mat[i,6], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),7] <- rep(lm_results_bothtime_mat[i,7], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),8] <- rep(lm_results_bothtime_mat[i,8], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),9] <- rep(lm_results_bothtime_mat[i,9], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),10] <- rep(lm_results_bothtime_mat[i,10], times = 20)
  model_response[as.numeric(k):as.numeric(k+19),11] <- rep(lm_results_bothtime_mat[i,11], times = 20)
  resp <- c(as.numeric(model_response[k,4])*VR_seq + as.numeric(model_response[k,7]))# calculate response line as y = mx + b
  model_response[as.numeric(k):as.numeric(k+19),12] <- resp
  
  k = k + 20
  
}
model_response <- as.data.frame(model_response)
model_response$predicted <- as.numeric(model_response$predicted)
model_response$Disturbance <- as.factor(model_response$Disturbance)
model_response$transience <- as.factor(model_response$transience)
model_response$Nitrogen <- as.factor(model_response$Nitrogen)
model_response$Slope <- as.numeric(model_response$Slope)
model_response$CI_lwr_slope <- as.numeric(model_response$CI_lwr_slope)
model_response$CI_upper_slope <- as.numeric(model_response$CI_upper_slope)
model_response$Intercept <- as.numeric(model_response$Intercept)
model_response$CI_lwr_int <- as.numeric(model_response$CI_lwr_int)
model_response$CI_upper_int <- as.numeric(model_response$CI_upper_int)

# nitrogen order
levels(model_response$Nitrogen)
model_response$Nitrogen <- factor(model_response$Nitrogen, levels = 
                                    c("0","1","2","3.4","5.4","9.5","17","27.2"))

# transient order
levels(model_response$transience)
model_response$transience <- factor(model_response$transience, 
                                    levels = c("Transient", "Post-transient"))

model_response$VR <- NA
model_response$VR <- rep(VR_seq, times = 32)

# produce figure
labels <- c("0" = "Intact in 1982", "1"= "Disturbed in 1982")

fig4_bothtime <- ggplot(data = vr_st_df_bothtime, aes(x=VR,y=Stability, col=Nitrogen))+
  geom_point(shape = 21, alpha = 0.4)+
  facet_grid(transience~ Disturbance, labeller=labeller(Disturbance = labels))+
  # geom_smooth(data = vr_st_df_bothtime, aes(x=VR, y=Stability,col=Nitrogen), method = 'lm', 
  #             se=FALSE, 
  #             fullrange = TRUE, 
  #             alpha = 0.8) +
  geom_line(data = model_response, aes(x = VR, y = predicted, col = Nitrogen),
            linewidth = 1.2) + # INSERT LINES FROM INTERCEPTS AND SLOPES IN THIS DATAFRAME
  geom_vline(aes(xintercept=1), linetype = "dashed")+
  labs(x="Synchrony", y="Stability")+
  theme_bw() +
  #labs(tag = "A") +
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
        text = element_text(size = 16))+
  ylim(0,8) + # should not be less than 0
  scale_colour_viridis_d(option = "D",direction=-1, na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))
fig4_bothtime

pdf(file = "Figures/Figure4_stabilityvssynch_bothtime.pdf",width = 8, height = 5)
fig4_bothtime
dev.off()

#### FigS2. slopes and intercepts ####
# For slopes
view(lm_results_bothtime)
lm_results_bothtime$Disturbance <- case_match(lm_results_bothtime$Disturbance, 
                                              "0" ~ "Intact in 1982", "1" ~ "Disturbed in 1982")
lm_results_bothtime$Disturbance <- factor(lm_results_bothtime$Disturbance,
                                          levels = c("Intact in 1982", "Disturbed in 1982"))

fig2_bothtime_slopes <- ggplot(data = lm_results_bothtime) + 
  geom_point(aes(x=Slope, y = Nitrogen, color = Nitrogen))+
  geom_errorbar(aes(xmin=CI_lwr_slope, xmax=CI_upper_slope, y = Nitrogen, color = Nitrogen)) +
  facet_grid(transience~ Disturbance)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        text = element_text(size = 16))+ # calculate baseline as vline (nitrogen == 0)
  labs(x = 'Slope', y = "Nitrogen addition", tag ="A") +
  scale_colour_viridis_d(option = "D",direction=-1, na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2")) +
  geom_vline(data = . %>% group_by(Disturbance) %>% distinct(), 
             aes(xintercept=baseline_slopes), color = 'black', linetype = 'dashed')
fig2_bothtime_slopes

# For intercepts

fig2_bothtime_int <- ggplot(data = lm_results_bothtime) + 
  geom_point(aes(x=Intercept, y = Nitrogen, color = Nitrogen))+
  geom_errorbar(aes(xmin=CI_lwr_int, xmax=CI_upper_int, y = Nitrogen, color = Nitrogen)) +
  facet_grid(transience~ Disturbance)+
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, face = "plain"),
        legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        text = element_text(size = 16))+
  labs(x = 'Intercept', y = "", tag = "B") +
  scale_colour_viridis_d(option = "D",direction=-1, na.value="grey72",
                         name="Nitrogen Addition",
                         breaks=c(0, 1, 2, 3.4, 5.4, 9.5, 17, 27.2),
                         labels=c("0.0", "1.0", "2.0", "3.4", "5.4", "9.5", "17.0", "27.2"))+
  geom_vline(data = . %>% group_by(Disturbance) %>% distinct(), 
             aes(xintercept=baseline_ints), color = 'black', linetype = 'dashed')
fig2_bothtime_int

# Full Figure 3 figure for collaborators
library(patchwork)
pdf(file = "Figures/SupFigure2_slopes_int.pdf",width = 10, height = 5)
(fig2_bothtime_slopes | fig2_bothtime_int)
dev.off()
