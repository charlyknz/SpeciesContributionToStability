## R Script to analyse AUC data ##
# Note: Please run 01SITES_createTidyData.R first 

####packages ####
library(tidyverse)
library(readxl)
library(cowplot)
library(here)

#### data ####
data <- zoo.stab.auc.prefin 
 # or 
data <- read.csv2('~/Desktop/phD/SpeciesContributionToStability/SITES_Data/AUCdata_2.csv', sep = ';')
str(data)

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

plot_grid( AUC.RR_mean.pi, AUC.pi_mean.pi,ncol = 2,hjust = -1,labels = c('(a)','(b)', '(c)', '(d)'),rel_heights = c(2,3))
ggsave(plot = last_plot(), width = 12, height = 5, file = here('output/Fig5new.png'))

#### correlation structure ####
str(Raw.dom.zoo)
Raw.dom.zoo$pos.AUC.RR <- abs(Raw.dom.zoo$mean.AUC.RR)
Raw.dom.zoo$pos.AUC.pi <- abs(Raw.dom.zoo$mean.AUC.pi)

#press 
corrPlot_press <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Press'), y = 'mean.AUC.RR', x = 'relat.dom',ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_press

corrPlot_pressPi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Press'), y = 'mean.AUC.pi', x = 'relat.dom',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pressPi

#pulse
corrPlot_pulse <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse'), y = 'mean.AUC.RR', x = 'relat.dom',ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulse

corrPlot_pulsePi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse'), y = 'mean.AUC.pi', x = 'relat.dom',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsePi

#pulsepress
corrPlot_pulsepress <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse & Press'), y = 'mean.AUC.RR', x = 'relat.dom',ylab = 'Absolute Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsepress

corrPlot_pulsepressPi <- ggscatter(subset(Raw.dom.zoo, Treatment == 'Pulse & Press'), y = 'mean.AUC.pi', x = 'relat.dom',ylab = 'Relative Contribution to stability', xlab = 'Relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot_pulsepressPi

SITES_corr <- plot_grid(corrPlot_press,corrPlot_pressPi,corrPlot_pulse,corrPlot_pulsePi,corrPlot_pulsepress,corrPlot_pulsepressPi, ncol = 2)
SITES_corr

ggsave(plot = SITES_corr, file = here('output/SITES_correlation.png'), width = 8, height = 12)

#### size ####

zoo.size <- read_csv("SITES_Data/zooplanktonMeanLength.csv") %>% 
  select(-...1  )

names(zoo.size)

## adjust Treatment labels
zoo.size$Treatment[zoo.size$Treatment == "F" ]<- "Pulse"
zoo.size$Treatment[zoo.size$Treatment == "FS"] <- "Pulse & Press"
zoo.size$Treatment[zoo.size$Treatment == "S"] <- "Press"

### join size with relat.dom data ###
zoo.with.size <- left_join(Raw.dom.zoo, zoo.size, by = c("Taxa", "Treatment" )) #%>%
names(zoo.with.size)
unique(zoo.with.size$Taxa)
hist(log(zoo.with.size$mean.size))

corrPlot <- ggscatter(zoo.with.size, y = 'mean.size', x = 'relat.dom',ylab = 'Mean size (in um)', xlab = 'Mean relative dominance',add = 'reg.line', conf.int = T)  +
  stat_cor( label.x = 0.1) 
corrPlot

### AUC plot SIZE ###

zoo.with.traits <- zoo.stab.auc.prefin %>%                                              # ranked by mean.pi
  gather('AUC.pi','AUC.RR', key = "AUC.stability", value = "AUC.value") %>%
  dplyr::group_by( Taxa, Treatment, AUC.stability) %>% 
  dplyr::mutate(mean.value = mean(AUC.value),
                sd.AUC = sd(AUC.value),
                se.value = sd.AUC/sqrt(n())) %>%
  mutate(trend = ifelse(AUC.value < 0, 'negative', 'positive')) %>%
  mutate(mean.trend = ifelse(mean.value < -0.1, 'negative', ifelse(mean.value > 0.1, 'positive', 'neutral'))) %>%
  left_join(., zoo.size, by = c("Taxa", "Treatment", 'Zooplankton_group' )) #%>%
names(zoo.with.traits)
unique(zoo.with.traits$Taxa)
hist(log(zoo.size$mean.size))

AUC.pi.size <- ggplot(subset(zoo.with.traits, AUC.stability == 'AUC.pi'), aes(x = mean.size, y=AUC.value, color = trend))+
  geom_point(alpha = .2)+
  geom_point(aes( x= mean.size, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = mean.size, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .1, color = 'black')+
  # geom_text(aes(x = mean.size,y = mean.value, label = Taxa), nudge_x = 0.05,nudge_y = 15,vjust = 0,color = "grey30" ,size = 2.5, angle = 85) +
  scale_color_manual(values = c('#D56060', '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  scale_x_log10(limits = c(40, 3200),breaks = c(100, 300, 1000, 3000))+
  geom_hline(yintercept = 0, alpha = 0.5) +
  labs(x = "Size (in µm on Ln-transformed scale)", y = 'Relative contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))#+
# theme(strip.text.x = element_text(color = 'white'),
#      strip.background =element_rect(fill = 'white', linetype = 0))
AUC.pi.size 

AUC.rr.size <- ggplot(subset(zoo.with.traits, AUC.stability == 'AUC.RR'), aes(x = mean.size, y=AUC.value, color = trend))+
  geom_point(alpha = .2)+
  geom_point(aes( x= mean.size, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = mean.size, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .1, color = 'black')+
  # geom_text(aes(x = mean.size,y = mean.value, label = Taxa), nudge_x = 0.1,nudge_y = 15,vjust = 0,color = "grey30" ,size = 2.5, angle = 85) +
  scale_color_manual(values = c('#D56060', '#68789E'))+
  scale_fill_manual(values = c('#D56060','grey', '#68789E'))+
  geom_hline(yintercept = 0, alpha = 0.5) +
  scale_x_log10()+
  labs(x = "Size (in µm on Ln-transformed scale)", y = 'Absolute contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "cm"))#+
# theme(strip.text.x = element_text(color = 'white'),
#        strip.background =element_rect(fill = 'white', linetype = 0))
AUC.rr.size

plot_grid(AUC.rr.size, corrPlot,AUC.pi.size,ncol = 2,hjust = -0.2,labels = c('(a)','(b)', '(c)', '(d)'),rel_widths = c(2.5,1.5))
ggsave(plot = last_plot(), width = 10, height = 7, file = here('output/Fig5Supplement.png'))


### Dominance Ranking ###

# 1.Rank
zoo.dummy <- data %>%
  group_by(Taxa) %>%
  summarise(mean.relat.dom = mean(mean.con.pi, na.rm = T)) %>%
  arrange(desc(mean.relat.dom))%>%
  mutate(ranking.by.mean.pi =  rank(dplyr::desc(mean.relat.dom))) %>%
  dplyr::arrange( ranking.by.mean.pi) 

# 2. Merge Rank with species contributions

zoo.effect <- left_join(data, zoo.dummy, by = c('Taxa')) %>%
  gather('AUC.pi','AUC.RR', key = "AUC.stability", value = "AUC.value") %>%
  dplyr::group_by( Taxa, Treatment, ranking.by.mean.pi, AUC.stability) %>% 
  dplyr::mutate(mean.value = mean(AUC.value),
                sd.AUC = sd(AUC.value),
                se.value = sd.AUC/sqrt(n())) %>%
  mutate(trend = ifelse(AUC.value < 0, 'negative', 'positive')) %>%
  mutate(mean.trend = ifelse(mean.value < -0.1, 'negative', ifelse(mean.value > 0.1, 'positive', 'neutral')))

str(zoo.effect)

### Plot ###
AUC.RR_dom <- ggplot(subset(zoo.effect, AUC.stability == 'AUC.RR'), aes(x = ranking.by.mean.pi, y= AUC.value, color = trend))+
  geom_hline(yintercept = 0, col="grey") +
  geom_point(alpha = .2)+
  geom_point(aes( x= ranking.by.mean.pi, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = ranking.by.mean.pi, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .8, color = 'black')+
  geom_text(aes(x = ranking.by.mean.pi,y = mean.value, label = Taxa), nudge_y = 20, vjust = 1,color = "grey30" ,size = 3, angle = 85) +
  scale_color_manual(values = c('#D56060',  '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))+
    labs(x = "Dominance ranking", y = 'Absolute contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(axis.title.x = element_text(size = 15,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.RR_dom 

AUC.pi_dom <- ggplot(subset(zoo.effect, AUC.stability == 'AUC.pi'), aes(x = ranking.by.mean.pi, y= AUC.value, color = trend))+
  geom_hline(yintercept = 0, col="grey") +
  geom_point(alpha = .2)+
  geom_point(aes( x= ranking.by.mean.pi, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = ranking.by.mean.pi, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .8, color = 'black')+
  geom_text(aes(x = ranking.by.mean.pi,y = mean.value, label = Taxa),  nudge_y = 10,vjust = 1,color = "grey30" ,size = 3, angle = 85) +
  scale_color_manual(values = c('#D56060', '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))+
   labs(x = "Dominance ranking", y = 'Relative contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(axis.title.x = element_text(size = 15,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.pi_dom 

plot_grid( AUC.RR_mean.pi, AUC.pi_mean.pi,AUC.RR_dom, AUC.pi_dom,ncol = 2,hjust = -1,labels = c('(a)','(b)', '(c)', '(d)'),rel_heights = c(2,3))
ggsave(plot = last_plot(), width = 14, height = 8, file = here('output/Fig5new.png'))


## look at relative dominance more closely ##

hist(zoo.effect$mean.relat.dom)
hist(Raw.dom.zoo$relat.dom)

# Ranking ~ mean of con.pi per Treatment 
MEAN.plot <- ggplot(zoo.effect, aes (x = ranking.by.mean.pi, y = mean.relat.dom))+
  geom_point()+
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))

raw.plot <- Raw.dom.zoo  %>%
  left_join(., zoo.dummy, by = c('Taxa')) %>%
  ggplot(., aes (x = ranking.by.mean.pi, y = relat.dom))+
  geom_point()+
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))
  
raw.plot.col <- Raw.dom.zoo  %>%
  left_join(., zoo.dummy, by = c('Taxa')) %>%
  ggplot(., aes (x = ranking.by.mean.pi, y = relat.dom, color = Lake, shape = Experiment))+
  geom_point()+
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))

blank <- ggplot()+
  geom_blank()+
  theme( panel.background = element_rect(fill = 'white', color = 'white'))
plot_grid(raw.plot, MEAN.plot, raw.plot.col, blank, rel_heights = c(2, 1.7))
ggsave(plot = last_plot(), file = here('output/RankingRelatDom.png'),width = 10, height = 8)


