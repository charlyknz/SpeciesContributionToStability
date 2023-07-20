## R Script to analyse AUC data ##
# Note: Please run 01SITES_createTidyData.R first 

####packages ####
library(tidyverse)
library(readxl)
library(cowplot)
library(here)
library(ggpubr)

#### data ####
data <- zoo.stab.auc.prefin # run 01 before
 # or 
data <- read.csv2('~/Desktop/phD/SpeciesContributionToStability/complete/AUCdata_3.csv', sep = ';')
str(data)

unique(data$Taxa)
### dominance and contribution to (in-)stability ###

#1. calculate Mean relative dominance & contributions for each Taxa in each Lake, Experiment
Raw.dom.zoo <- data %>%
  group_by(Lake, Experiment, Taxa) %>%
  mutate(relat.dom = mean(mean.con.pi, na.rm = T)) %>%
  group_by(Lake, Experiment, Taxa, relat.dom, Treatment) %>%
  summarise(mean.AUC.RR = mean(AUC.RR),
            sd.AUC.RR = sd(AUC.RR),
            se.AUC.RR = sd.AUC.RR/sqrt(n()),
            mean.AUC.pi = mean(AUC.pi),
            sd.AUC.pi = sd(AUC.pi),
            se.AUC.pi = sd.AUC.pi/sqrt(n())) %>%
  mutate(mean.trend.RR = ifelse(mean.AUC.RR < -0.1, 'negative', ifelse(mean.AUC.RR > 0.1, 'positive', 'neutral'))) %>%
  mutate(mean.trend.pi = ifelse(mean.AUC.pi < -0.1, 'negative', ifelse(mean.AUC.pi > 0.1, 'positive', 'neutral'))) 

# plot output
AUC.RR_mean.pi <- ggplot(Raw.dom.zoo, aes(x = relat.dom, y= mean.AUC.RR, color = mean.trend.RR,  shape = Experiment))+
  geom_hline(yintercept = 0, col="grey") +
  geom_point(alpha = .6, size = 2.5)+
  scale_color_manual(values = c('#D56060', 'grey','#68789E'))+
 # scale_x_reverse()+
  scale_y_continuous(limits = c(-22,22), breaks = seq(-20,20,10))+
  labs(x = "Relative dominance", y = 'Absolute contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "cm"))+
  theme(axis.title.x = element_text(size = 15,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.RR_mean.pi 

AUC.pi_mean.pi <- ggplot(Raw.dom.zoo, aes(x = relat.dom, y= mean.AUC.pi, color = mean.trend.pi, shape = Experiment))+
  geom_hline(yintercept = 0, col="grey") +
  geom_point(alpha = .6, size = 2.5)+
  scale_color_manual(values = c('#D56060', 'grey','#68789E'))+
  scale_y_continuous(limits = c(-12,21), breaks = c(-10,0,10,20))+
  # scale_x_reverse()+
  labs(x = "Relative dominance", y = 'Relative contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "cm"))+
  theme(axis.title.x = element_text(size = 15,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.pi_mean.pi 

cowplot::plot_grid( AUC.RR_mean.pi, AUC.pi_mean.pi,ncol = 2,hjust = -1,labels = c('(a)','(b)', '(c)', '(d)'),rel_heights = c(2,3))
ggsave(plot = last_plot(), width = 12, height = 5, file = here('complete/Fig5new.png'))



#### correlation structure ####
str(Raw.dom.zoo)
Raw.dom.zoo$pos.AUC.RR <- abs(Raw.dom.zoo$mean.AUC.RR)
Raw.dom.zoo$pos.AUC.pi <- abs(Raw.dom.zoo$mean.AUC.pi)

#press 
corrPlot_press <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Press'), y = 'pos.AUC.RR', x = 'relat.dom',cor.method = 'spearman',
                            ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_press

corrPlot_pressPi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Press'), y = 'pos.AUC.pi', x = 'relat.dom',cor.method = 'spearman',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pressPi

#pulse
corrPlot_pulse <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse'), y = 'pos.AUC.RR', x = 'relat.dom',cor.method = 'spearman',ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulse

corrPlot_pulsePi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse'), y = 'pos.AUC.pi', x = 'relat.dom',cor.method = 'spearman',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsePi

#pulsepress
corrPlot_pulsepress <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse & Press'), y = 'pos.AUC.RR', x = 'relat.dom',cor.method = 'spearman',ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsepress

corrPlot_pulsepressPi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse & Press'), y = 'pos.AUC.pi', x = 'relat.dom',cor.method = 'spearman',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsepressPi

SITES_corr <- cowplot::plot_grid(corrPlot_press,corrPlot_pressPi,corrPlot_pulse,corrPlot_pulsePi,corrPlot_pulsepress,corrPlot_pulsepressPi, ncol = 2)
SITES_corr

ggsave(plot = SITES_corr, file = here('complete/SITES_correlation.png'), width = 8, height = 12)
