## Script to analyse simulated disturbance effects on community ## 
# by Charlotte Kunze

####load packages ####
#library(reshape2)
library(tidyverse)
library(ggpubr)
library(psych)
library(cowplot)
library(here)


# import data
Stab_Alpha_AUC_M5<- read.csv('BEFD_createdData/StabAlphaAUC.csv')
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

rectLim1 <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                    x1 =c(-4,NA, NA) ,
                    x2 = c(-0.15,NA, NA) ,
                    y1 = c(3.5,NA, NA) ,
                    y2 = c(5,NA, NA) )


TextLim1  <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                        x = -3.98,
                        y = 4,
                        text = c('All equally sensitive','',''))


Lim1RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit1'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  geom_point(alpha= .3)+
  geom_text(data = TextLim1, inherit.aes = F, aes(x = x, y = y, label = text), size = 3.4, hjust = 0) +
  geom_rect(data = rectLim1, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
  scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
  scale_y_continuous(limits = c(-6, 6), breaks = c(-4,-2,0,2,4))+
  scale_x_continuous(limits = c(-4, 4), breaks = c(-2,0,2))+
  labs(x = '',y = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.background =element_rect(),
        strip.text = element_text(size = 12, face = 'bold'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim1RRpi


rectLim2 <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                   x1 =c(-4,NA, NA) ,
                   x2 = c(-0.15,NA, NA) ,
                   y1 = c(17,NA, NA) ,
                   y2 = c(23,NA, NA) )


TextLim2  <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                    x = -3.98,
                    y = 20,
                    text = c('One tolerant sp.','',''))


Lim2RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit2'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
 scale_y_continuous(limits = c(-27, 27), breaks = c(-20,-10,0,10,20))+
  scale_x_continuous(limits = c(-4, 4), breaks = c(-2,0,2))+
  geom_point(alpha= .3)+
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  geom_text(data = TextLim2, inherit.aes = F, aes(x = x, y = y, label = text), size = 3.5, hjust = 0) +
  geom_rect(data = rectLim2, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
  labs(x = '', y = "Absolute contribution to stability")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  theme(legend.position="none")
Lim2RRpi


rectLim3 <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                   x1 =c(-4,NA, NA) ,
                   x2 = c(-0.15,NA, NA) ,
                   y1 = c(26,NA, NA) ,
                   y2 = c(32.99,NA, NA) )


TextLim3  <- tibble(Model = c('Press', 'Pulse', 'Pulse & Press'),
                    x = -3.98,
                    y = 30,
                    text = c('One sensitive sp.','',''))

Lim3RRpi<-ggplot(subset(stab.alpha, Limit == 'Limit3'), aes(x=AUC.pi, y=AUC.RR,
                                                                  col=species,
                                                                  size= relAlpha)) +
  geom_point(alpha= .3)+
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  geom_text(data = TextLim3, inherit.aes = F, aes(x = x, y = y, label = text), size = 3.5, hjust = 0) +
  geom_rect(data = rectLim3, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
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

par(mar=c (5.1, 4.1, 4.1, 2.1))
ggarrange(Lim1RRpi, Lim2RRpi, Lim3RRpi,ncol = 1, nrow = 3, vjust = 3.0,hjust=-1.5,widths = c(1,1), labels = c( '(a)', '(b)', '(c)'))
ggsave(plot = last_plot(), width = 9, height =10, file = here('output/Fig.2.pdf'))

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
  labs(y='Absolute contribution to stability', x = "Relative competitiveness")+
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
  labs(y= 'Relative contribution to stability', x = "Relative competitiveness")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(strip.text.x = element_text(color = 'white'),
        strip.background =element_rect(fill = 'white', linetype = 0))+
  theme(text = element_text(size = 16)) +  theme(legend.position="")

Lim3pi#ggsave(plot = RRpi, 'RRPi.png',width = 8, height = 4)

par(mar=c (5.1, 4.1, 4.1, 2.1))
ggarrange(Lim1R,Lim1pi,  Lim2R,Lim2pi,Lim3R,Lim3pi,hjust = -1, ncol = 2, nrow = 3, labels =  c( '(a)',' ', '(b)',' ', '(c)'))
ggsave(plot = last_plot(), width = 12, height = 12, file = here('output/Supplement_AllDist.png'))


