#### RAD analysis ####
# Charlotte Kunze 20.01.22

# NOTE: you need to run 01AUC_Zoo_PabloData.R first to create the data needed here!!


#### packages ####
library(tidyverse)
library(nlme)
library(lmerTest)
library(performance)
library(see)
library(readxl)
library("ggsci")
library(ggpubr)

#### check data ####
names(zoo.stab.auc.prefin)
str(zoo.stab.auc.prefin)
names(zoo.data3)

# check data levels
levels(as.factor(zoo.stab.auc.prefin$Taxa))
sapply(lapply(zoo.stab.auc.prefin, unique), length)

hist(zoo.stab.auc.prefin$AUC.pi)
hist(zoo.stab.auc.prefin$AUC.RR)

#change variables to factors for analysis
zoo.stab.auc.prefin$Experiment <- as.factor(zoo.stab.auc.prefin$Experiment)
zoo.stab.auc.prefin$Treatment <- as.factor(zoo.stab.auc.prefin$Treatment)
zoo.stab.auc.prefin$Lake <- as.factor(zoo.stab.auc.prefin$Lake)


#remove entries with NA
which(is.na(zoo.stab.auc.prefin))
zoo.stab.auc.prefin  <-zoo.stab.auc.prefin%>%
  drop_na(AUC.RR, AUC.pi)


#### Rank Abundance Diagrams ####

### barplot of species' mean.pi pre- and past-disturbance ###
zoo.data3 %>% 
  select(Lake, Experiment, Exp_day, Treatment, Taxa, mean.pi) %>%
  mutate(disturbance = paste(ifelse(Exp_day < 7, 'pre-dist', 'past-dist'))) %>%
  dplyr::group_by(Taxa, disturbance) %>% 
  dplyr::summarise(mean = mean(mean.pi))%>%
  arrange(desc(mean))%>%
  filter(mean >0.1) %>%
  ggplot(., aes(x = disturbance, y = mean, fill = Taxa))+
  geom_col(color = 'black')
### ###


### RANK after most dominant species in control over time (con.pi) ### 

# 1. ranking
zoo.dummy.con <- zoo.data3 %>% 
  select(Lake, Experiment, Exp_day, Treatment, Taxa, con.pi, Abundance, Biomass, USI) %>%
  mutate(one = as.numeric(paste(1)) )%>%
  add_count(USI, wt = one) %>% # count number of occurrences in data
  dplyr::group_by(Taxa) %>% 
  filter(n > 2) %>%  #remove species which are not present in stab.auc dataset
  filter(!Taxa %in% c("Nauplii", "Small cladoceran")) %>%
  dplyr::summarise(mean = mean(con.pi))%>%
  mutate(ranking.by.mean.pi =  rank(dplyr::desc(mean))) %>%
  dplyr::arrange( ranking.by.mean.pi) 

# 2. join datasets and calculate mean auc
zoo.effect.plot <- left_join(zoo.stab.auc.prefin, zoo.dummy.con, by = c("Taxa"))  %>%                                              # ranked by mean.pi
  gather('AUC.pi','AUC.RR', key = "AUC.stability", value = "AUC.value") %>%
  dplyr::group_by(Taxa, Treatment, ranking.by.mean.pi, AUC.stability) %>% 
  dplyr::mutate(mean.value = mean(AUC.value),
                sd.AUC = sd(AUC.value),
                se.value = sd.AUC/sqrt(n())) %>%
  mutate(trend = ifelse(AUC.value < 0, 'negative', 'positive')) %>%
  mutate(mean.trend = ifelse(mean.value < -0.1, 'negative', ifelse(mean.value > 0.1, 'positive', 'neutral')))

names(zoo.effect.plot)

#### zoo traits ####
zoo.feeding <- read_excel("SITES_Data/collectedZooplanktonFeedingMode.xlsx")
names(zoo.feeding)

zoo.feeding$Treatment[zoo.feeding$Treatment == "F" ]<- "Pulse"
zoo.feeding$Treatment[zoo.feeding$Treatment == "FS"] <- "Pulse & Press"
zoo.feeding$Treatment[zoo.feeding$Treatment == "S"] <- "Press"

unique(zoo.feeding$consumption_mode)
levels(as.factor(zoo.effect.plot$mean.trend))

#### size data ####
zoo.size <- read_csv("SITES_Data/zooplanktonMeanLength.csv") %>% 
  select(-...1  )
names(zoo.size)
zoo.size$Treatment[zoo.size$Treatment == "F" ]<- "Pulse"
zoo.size$Treatment[zoo.size$Treatment == "FS"] <- "Pulse & Press"
zoo.size$Treatment[zoo.size$Treatment == "S"] <- "Press"

#### join size with RAD data ####
names(zoo.effect.plot)
names(zoo.feeding)
zoo.with.traits <- left_join(zoo.effect.plot, zoo.feeding, by =c("Taxa", "Treatment",'Zooplankton_group' ))
zoo.with.traits <- left_join(zoo.with.traits, zoo.size, by = c("Taxa", "Treatment", 'Zooplankton_group' )) #%>%

names(zoo.with.traits)

#### correlation plot ####
corr_size_AUCpi <- zoo.with.traits %>%
  mutate(log.size = log10(mean.size)) %>%
  filter(AUC.stability != 'AUC.RR') %>%
ggscatter(., x = 'log.size', y = 'AUC.value',add = 'reg.line', conf.int = T, xlab = 'size (logartihmized scale)', ylab='Relative contribution to stability' ) +
  stat_cor( label.x = 1.8)           # Add correlation coefficient

corr_size_AUCrr <- zoo.with.traits %>%
  mutate(log.size = log10(mean.size)) %>%
  filter(AUC.stability == 'AUC.RR') %>%
  ggscatter(., x = 'log.size', y = 'AUC.value',add = 'reg.line', conf.int = T, xlab = 'size (logartihmized scale)', ylab='Absolute contribution to stability' ) +
  stat_cor( label.x = 1.8)           # Add correlation coefficient

plot_grid(corr_size_AUCpi, corr_size_AUCrr, labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = 'Supplement_correlation_size_stability.png', width = 8, height = 4)

#### Dominance plot size ####

### Fig. 5 a,b ###
AUC.RR_dom <- ggplot(subset(zoo.with.traits, AUC.stability == 'AUC.RR'), aes(x = ranking.by.mean.pi, y= AUC.value, color = trend))+
  geom_point(alpha = .2)+
  geom_point(aes( x= ranking.by.mean.pi, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = ranking.by.mean.pi, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .8, color = 'black')+
  geom_text(aes(x = ranking.by.mean.pi,y = mean.value, label = Taxa), nudge_y = 20, vjust = 1,color = "grey30" ,size = 3, angle = 85) +
  scale_color_manual(values = c('#D56060',  '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  geom_hline(yintercept = 0, alpha = 0.8) +
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))+
  # scale_color_manual(values = c('#A9A9A9', '#000000'))+
  # scale_shape_manual(values = c(16, 1))+
  labs(x = "Ranking after relat. dominance", y = 'Absolute contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.RR_dom 

AUC.pi_dom <- ggplot(subset(zoo.with.traits, AUC.stability == 'AUC.pi'), aes(x = ranking.by.mean.pi, y= AUC.value, color = trend))+
  geom_point(alpha = .2)+
  geom_point(aes( x= ranking.by.mean.pi, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = ranking.by.mean.pi, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .8, color = 'black')+
  geom_text(aes(x = ranking.by.mean.pi,y = mean.value, label = Taxa),  nudge_y = 10,vjust = 1,color = "grey30" ,size = 3, angle = 85) +
  scale_color_manual(values = c('#D56060', '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  geom_hline(yintercept = 0, alpha = .8) +
  scale_x_continuous(limits = c(0.5,18.5 ), breaks = seq(1,18,2))+
  # scale_color_manual(values = c('#A9A9A9', '#000000'))+
  # scale_shape_manual(values = c(16, 1))+
  labs(x = "Ranking after relat. dominance", y = 'Relative contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 16, face = 'bold'))
AUC.pi_dom 

### Fig. 5 c,d ###
AUC.pi.size <- ggplot(subset(zoo.with.traits, AUC.stability == 'AUC.pi'), aes(x = mean.size, y=AUC.value, color = trend))+
  geom_point(alpha = .2)+
  geom_point(aes( x= mean.size, y= mean.value, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = mean.size, ymin = mean.value-se.value, ymax = mean.value+se.value,color = mean.trend), width = .1, color = 'black')+
 # geom_text(aes(x = mean.size,y = mean.value, label = Taxa), nudge_x = 0.05,nudge_y = 15,vjust = 0,color = "grey30" ,size = 2.5, angle = 85) +
  scale_color_manual(values = c('#D56060', '#68789E'))+
  scale_fill_manual(values = c('#D56060', 'grey','#68789E'))+
  scale_x_log10(limits = c(40, 3200),breaks = c(100, 300, 1000, 3000))+
  geom_hline(yintercept = 0, alpha = 0.5) +
  labs(x = "size (in µm on log-transformed scale)", y = 'Relative contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))
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
  labs(x = "size (in µm on log-transformed scale)", y = 'Absolute contribution to stability') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "cm"))+
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))
AUC.rr.size

### all Fig. 5 ###
plot_grid(AUC.RR_dom, AUC.pi_dom, AUC.rr.size, AUC.pi.size,ncol = 2,hjust = -1.1,labels = c('(a)','(b)', '(c)', '(d)'),rel_widths = c(1,1))
ggsave(plot = last_plot(), width = 12, height = 8, file = 'Fig5.png')


### Supplement Fig. 5 ###
levels(zoo.with.traits$consumption_mode) <- c("omnivore", "Omnivore, carnivore" ,"Omnivore, herbivore",  "herbivore","zooplanktivore")

AUC.rr.ft <- zoo.with.traits %>%
  filter( AUC.stability == 'AUC.RR') %>%
  group_by(consumption_mode, AUC.stability) %>%
  mutate(mean = mean(AUC.value, na.rm  = T), 
          sd = sd(AUC.value, na.rm  = T),
         se = sd/sqrt(n())) %>%
  mutate(mean.trend = ifelse(mean < 0, 'negative', ifelse(mean> 0.1,'positive', 'neutral')))%>%
  ggplot(., aes(x = consumption_mode, y=AUC.value, color = trend))+
  geom_point(alpha = .3)+
  geom_point(aes( x= consumption_mode, y= mean, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = consumption_mode, ymin = mean-se, ymax = mean+se,color = mean.trend), width = .2, color = 'black')+
  scale_color_manual(values = c('#D56060', '#68789E', 'grey'))+
  scale_fill_manual(values = c('#D56060', '#68789E', 'grey'))+
  geom_hline(yintercept = 0, alpha = 0.5) +
  labs(x = "Feeding mode", y = 'AUC.RR (Functional Stability)') +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 13, face = 'bold'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

AUC.rr.ft

AUC.pi.ft<- zoo.with.traits %>%
  filter( AUC.stability == 'AUC.pi') %>%
  group_by(consumption_mode, AUC.stability) %>%
  mutate(mean = mean(AUC.value, na.rm  = T), 
         sd = sd(AUC.value, na.rm  = T),
         se = sd/sqrt(n())) %>%
  mutate(mean.trend = ifelse(mean < -0.1, 'negative', ifelse(mean> 0.1,'positive', 'neutral')))%>%
  ggplot(., aes(x = consumption_mode, y=AUC.value, color = trend))+
  geom_point(alpha = .3)+
  geom_point(aes( x= consumption_mode, y= mean, fill = mean.trend), size = 3, pch = 21, color = 'black') +
  geom_errorbar(aes(x = consumption_mode, ymin = mean-se, ymax = mean+se,color = mean.trend), width = .2, color = 'black')+
  scale_color_manual(values = c('grey', '#68789E'))+
  scale_fill_manual(values = c('grey','#68789E'))+
  geom_hline(yintercept = 0, alpha = 0.5) +
  labs(x = "Feeding mode", y =  expression(AUC.~Delta~'pi'~'('~Compositional~Stability~')')) +
  theme_bw() +
  facet_grid(~Treatment)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        #legend.background = element_rect(linetype = "solid", colour = "grey", size = 1),
        legend.key.size = unit(1, "cm"))+
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 13, face = 'bold'))+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

AUC.pi.ft
plot_grid(AUC.rr.ft, AUC.pi.ft, ncol = 2,hjust = -1.1,labels = c('(a)','(b)', '(c)', '(d)'),rel_widths = c(1,1))
ggsave(plot = last_plot(), width = 9, height = 4, file = 'Supplement_AUCfeedingMode.png')


