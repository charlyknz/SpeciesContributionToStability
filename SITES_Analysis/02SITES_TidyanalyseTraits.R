## R Script to analyse AUC data ##
# Note: Please run 01SITES_createTidyData.R first 

####packages ####
library(tidyverse)
library(readxl)
library(cowplot)

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
  scale_y_continuous(limits = c(-30,30), breaks = seq(-30,30,10))+
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
  scale_y_continuous(limits = c(-13,23), breaks = c(-10,0,10,20))+
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
ggsave(plot = last_plot(), width = 14, height = 8, file = here('output/Fig5new.png'))


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
