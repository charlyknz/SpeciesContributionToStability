## Script to analyse simulated disturbance effects on community ## 
# by Charlotte Kunze

####load packages ####
#library(reshape2)
library(tidyverse)
library(nlme)
library(ggpubr)
library(psych)
library(cowplot)


# set working directory 
setwd("~/Desktop/BEFD22/BEFDisturbance/Untitled/modelOutput")


# import data
Stab_Alpha_AUC_M5 <- max
Stab_Alpha_AUC_M5<- read.csv('StabAlphaAUC.csv')
str(Stab_Alpha_AUC_M5)

#### Data wrangling ####
stab.alpha <- Stab_Alpha_AUC_M5 %>%
  mutate(maximum = max(competition), #calculate maximum alpha 
        relAlpha= maximum/competition) %>%# to measure relative competitiveness calculate 1/relAlpha
  distinct(species,relAlpha,Limit, sensitivity, runNumber, Model, AUC.RR, AUC.pi) %>%
  mutate(runID = paste(Model, runNumber, sep = '_')) 
  
names(stab.alpha)

stab.alpha$Model[stab.alpha$Model == 'press']<- 'Press'
stab.alpha$Model[stab.alpha$Model == 'pulse']<- 'Pulse'
stab.alpha$Model[stab.alpha$Model == 'pulsepress']<- 'Pulse & Press'

#### plots ####

# change in biomass and alpha
Lim1RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit1'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  geom_point(alpha= .3)+
  scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
  scale_y_continuous(limits = c(-6, 6), breaks = c(-4,-2,0,2,4))+
  scale_x_continuous(limits = c(-3, 3), breaks = c(-2,0,2))+
  labs(x = '',y = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.background =element_rect(),
        strip.text = element_text(size = 12, face = 'bold'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim1RRpi

unique(stab.alpha$species)
Lim2RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit2'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
 scale_y_continuous(limits = c(-27, 27), breaks = c(-20,-10,0,10,20))+
  scale_x_continuous(limits = c(-4, 4), breaks = c(-3,0,3))+
  geom_point(alpha= .3)+
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  labs(x = '', y = "Absolute contribution to stability")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim2RRpi

Lim3RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit3'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  geom_point(alpha= .3)+
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  scale_color_manual(values = c('#D56060', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_y_continuous(limits = c(-33, 33), breaks = c(-30,-15,0,15,30))+
  scale_x_continuous(limits = c(-4, 4), breaks = c(-2,0,2))+
  labs(x = 'Relative contribution to stability', y = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +  theme(legend.position="none")

Lim3RRpi
#ggsave(plot = RRpi, 'RRPi.png',width = 8, height = 4)

par(mar=c (5.1, 4.1, 4.1, 2.1))
ggarrange(Lim1RRpi, Lim2RRpi, Lim3RRpi,ncol = 1, nrow = 3, vjust = 3.7,hjust=-1.5,widths = c(1,1))
ggsave(plot = last_plot(), width = 9, height = 9, file = 'DistTypes_auc.png')

hist(stab.alpha$relAlpha)

#### Supplement Figure: relative competitiveness plot ####
Lim1R<-ggplot(subset(stab.alpha, Limit == 'Limit1'), aes(x=relAlpha, y=AUC.RR,
                                                            col=species)) +
  geom_hline(yintercept=0, col="grey")+
  geom_point(alpha= .3, size = 2)+
  scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
  labs(y = 'Absolute contribution to stability', x = "")+
  facet_wrap(~Model, ncol = 3)+
  scale_y_continuous(limits = c(-1.7, 1.1), breaks = c(-1.5,-1,-0.5,0,0.5,1))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  theme_bw() +
  theme(strip.background =element_rect(),
        strip.text = element_text(size = 12, face = 'bold'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim1R

Lim1pi<-ggplot(subset(stab.alpha, Limit == 'Limit1'), aes(x=relAlpha, y=AUC.pi,
                                                          col=species)) +
  geom_hline(yintercept=0, col="grey")+
  geom_point(alpha= .3, size = 2)+
  scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
  scale_y_continuous(limits = c(-0.1, 0.1), breaks = c(-0.1,-0.05,0,0.05,0.1))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  labs(y = 'Relative contribution to stability', x = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.background =element_rect(),
        strip.text = element_text(size = 12, face = 'bold'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim1pi

unique(stab.alpha$species)
Lim2R<-ggplot(subset(stab.alpha, Limit == 'Limit2'), aes(x=relAlpha, y=AUC.RR,
                                                            col=species)) +
  scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_y_continuous(limits = c(-33, 33), breaks = c(-30,-15,0,15,30))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  geom_point(alpha= .3, size = 2)+
  geom_hline(yintercept=0, col="grey")+
  labs(y = 'Absolute contribution to stability', x = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim2R

Lim2pi<-ggplot(subset(stab.alpha, Limit == 'Limit2'), aes(x=relAlpha, y=AUC.pi,
                                                         col=species)) +
  scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_y_continuous(limits = c(-4, 4), breaks = c(-4,-2,0,2,4))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  geom_point(alpha= .3, size = 2)+
  geom_hline(yintercept=0, col="grey")+
  labs(y = 'Relative contribution to stability', x = "")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim2pi

Lim3R<-ggplot(subset(stab.alpha, Limit == 'Limit3'), aes(x=relAlpha, y=AUC.RR,
                                                            col=species)) +
  geom_point(alpha= .3, size = 2)+
  geom_hline(yintercept=0, col="grey")+
  scale_color_manual(values = c('#D56060', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  scale_y_continuous(limits = c(-30, 25), breaks = c(-30,-20,-10,0,10,20))+
  labs(y='Absolute contribution to stability', x = "relative competitiveness")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(text = element_text(size = 16)) +  theme(legend.position="none")

Lim3R
Lim3pi<-ggplot(subset(stab.alpha, Limit == 'Limit3'), aes(x=relAlpha, y=AUC.pi,
                                                         col=species)) +
  geom_point(alpha= .3, size = 2)+
  geom_hline(yintercept=0, col="grey")+
  scale_color_manual(values = c('#D56060', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_x_continuous(limits = c(0.99, 1.17), breaks = c(1,1.05,1.10,1.15))+
  scale_y_continuous(limits = c(-3, 2.5), breaks = c(-3,-2,-1,0,1,2))+
  labs(y= 'Relative contribution to stability', x = "relative competitiveness")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(text = element_text(size = 16)) +  theme(legend.position="")

Lim3pi#ggsave(plot = RRpi, 'RRPi.png',width = 8, height = 4)

par(mar=c (5.1, 4.1, 4.1, 2.1))
ggarrange(Lim1R,Lim1pi,  Lim2R,Lim2pi,Lim3R,Lim3pi,hjust = -1, ncol = 2, nrow = 3)
ggsave(plot = last_plot(), width = 12, height = 12, file = 'AllDist_sil2.png')


