## R Script to analyse AUC data ##
# Note: Please run 01SITES_createTidyData.R first 

####packages ####
library(tidyverse)
library(cowplot)
library(here)
library(ggpubr)

#### data ####
data <- zoo.stab.auc.prefin #works only if you ran 01SITES_createTidyData.R
# or 
data <- read.csv2('~/Desktop/phD/SpeciesContributionToStability/SITES_Data/AUCdata_2.csv', sep = ';')
str(data)


#### Distance calculation ####

trans.data <- data %>%
  group_by(Lake, Experiment, Taxa) %>%
  mutate(relat.dom = mean(mean.con.pi, na.rm = T))%>%
  group_by(Lake, Experiment, Taxa, relat.dom, Treatment) %>%
  summarise(mean.AUC.RR = mean(AUC.RR), mean.AUC.pi = mean(AUC.pi))%>%
  mutate(mean.trend.RR = ifelse(mean.AUC.RR < -0.1, 'negative', ifelse(mean.AUC.RR > 0.1, 'positive', 'neutral'))) %>%
  mutate(mean.trend.pi = ifelse(mean.AUC.pi < -0.1, 'negative', ifelse(mean.AUC.pi > 0.1, 'positive', 'neutral')))%>%
   ungroup()%>%
  mutate(max.rr = max(abs(mean.AUC.RR)),
         max.pi = max(abs(mean.AUC.pi))) %>%
  ungroup() %>%
  mutate(AUC.trans.RR = mean.AUC.RR/max.rr,
         AUC.trans.pi = mean.AUC.pi/ max.pi)  %>%
  mutate(dist = sqrt( ( abs(AUC.trans.pi)^2+ abs(AUC.trans.RR)^2))) %>%
  drop_na(dist)  # 

# plot output
distance <- ggplot(trans.data, aes(x = relat.dom, y= dist, color = Taxa))+
  geom_hline(yintercept = 0, col="grey") +
  geom_point(alpha = .6, size = 2.5)+
  #scale_color_manual(values = c('#D56060', 'grey','#68789E'))+
  #scale_x_reverse()+
  #scale_y_continuous(limits = c(-22,22), breaks = seq(-20,20,10))+
  labs(x = "Relative dominance", y = 'Distance') +
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
distance 
ggsave(plot = last_plot(), file = here('output/DistanceRelDom.png'), width = 10, height = 5)

#### correlation ####
unique(trans.data$Treatment)
corrPlot_dist.press <- ggscatter(subset(trans.data, Treatment == 'Press'), y = 'dist', x = 'relat.dom',ylab = 'Distance', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_dist.press

corrPlot_dist.pulse <- ggscatter(subset(trans.data, Treatment == 'Pulse'), y = 'dist', x = 'relat.dom',ylab = 'Distance', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_dist.pulse

corrPlot_dist.presspulse <- ggscatter(subset(trans.data, Treatment == 'Pulse & Press'), y = 'dist', x = 'relat.dom',ylab = 'Distance', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_dist.presspulse

plot_grid(corrPlot_dist.press,corrPlot_dist.pulse,corrPlot_dist.presspulse, ncol = 1,hjust = -1.8,labels = c('press', 'pulse', 'pulsepress'))
ggsave(plot = last_plot(), file = here('output/correlation_distanceRelDomH.png'), width = 4, height = 12)


#### magnitude of cont ~ distance ####

corrPlot_pulse <- trans.data %>%
  filter(Treatment == 'Pulse') %>%
  mutate(magnitudeRR = abs(mean.AUC.RR)) %>%
  ggscatter(., y = 'dist', x = 'magnitudeRR',ylab = 'Distance', xlab = 'Magnitude Absolute Contribution',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulse

corrPlot_pulsePi <- trans.data %>%
  filter(Treatment == 'Pulse') %>%
  mutate(magnitudepi = abs(mean.AUC.pi)) %>%
  ggscatter(., y = 'dist', x = 'magnitudepi',ylab = 'Distance', xlab = 'Magnitude Relative Contribution',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsePi
ggsave(corrPlot_pulsePi, file = here('output/distanceRelContribution.png'), width = 4, height = 4)
