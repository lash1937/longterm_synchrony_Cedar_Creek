#####################################
# This script produces figure S?: comparing robustness of results for 7 year windows
# vs. 10 year windows
#####################################

# packages necessary
library(dplyr)
library(ggplot2)

# load in data
ten_PT <- readRDS(here::here("data/SEM_10yr_posttransient.rds"))
ten_T <- readRDS(here::here("data/SEM_10yr_transient.rds"))
seven_PT <- readRDS(here::here("data/SEM_posttransient.rds"))
seven_T <- readRDS(here::here("data/SEM_transient.rds"))

# combine all dataframes into one
ten_PT$window <- "10"
ten_PT$phase <- "Post-transient"

ten_T$window <- "10"
ten_T$phase <- "Transient"

seven_PT$window <- "7"
seven_PT$phase <- "Post-transient"

seven_T$window <- "7"
seven_T$phase <- "Transient"


window_comparison <- rbind(ten_PT, ten_T, seven_PT, seven_T)


# create new 'pathway' variable to compare windows across all pathways tested
window_comparison <- window_comparison %>%
  tidyr::unite("pathway", c("lhs", "rhs"), remove = FALSE, sep="_")

# define pathways of interest
var_tested <- c("TStability", "TVR", "TRichness", "TEvenness", "Nitrogen")

window_comparison <- window_comparison %>%
  filter(lhs %in% var_tested) %>%
  filter(rhs %in% var_tested) %>%
  filter(op == "~")
window_comparison$pathway <- as.character(window_comparison$pathway)
window_comparison$pathway <- as.factor(window_comparison$pathway)
window_comparison$window <- as.factor(window_comparison$window)
window_comparison$y <- 0


windows_test <- ggplot(data = window_comparison, aes(x = pathway, y = y)) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper, width = 0.2, color = window), position = position_dodge(width = 0.8)) +
  theme_minimal()+
  scale_color_manual(values = c("#D55E00", "skyblue"),
                    name="Window",
                    labels = c("10 years", "7 years"))+
  facet_wrap(~phase + group)+
  ylab("")+
  xlab("Pathway")+
  theme(axis.text.x = element_text(color = "grey20", size = 12,
                                   angle = 90, hjust = -0.12, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = .5, vjust = 0,
                                   face = "plain"),
        axis.title.x = element_text(color = "black", size = 12,
                                    angle = 0, hjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12,
                                    angle = 90, hjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12,
                                   angle = 0, hjust = 0, face = "plain"))




