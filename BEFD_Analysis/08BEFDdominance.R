#### R script for dominance plots ####
# by Charlotte Kunze 

#load packages and functions
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggpmisc) #for regression lines
library(here)

#### stabAlphaAUC dataset ####
max<- read.csv('BEFD_createdData/StabAlphaAUC.csv')

# data wrangling - remove duplicates created in the AUC loop
stab.auc <- distinct(max, Limit, Model, sensitivity, competition, relAlpha, max, runNumber, species, mean.treat.pi, mean.con.pi, AUC.pi, AUC.RR, AUC.totRR) %>%
  group_by(Limit, runNumber, Model) %>%
  mutate(MaxAlpha = max(competition)) %>%
  ungroup() %>%
  mutate(inv_relAlpha = MaxAlpha/competition) 

#compare measures
hist(stab.auc$mean.con.pi )
hist(stab.auc$inv_relAlpha)


#### correlations ####

## Relative Contributions ##
# press 
Lim1press <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit1') %>%
ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1press

Lim2press <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2press

Lim3press <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3press


# pulse 
Lim1pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1pulse

Lim2pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2pulse

Lim3pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3pulse


# pulsepress 
Lim1pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1pulsepress

Lim2pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2pulsepress

Lim3pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3pulsepress

corrRela<- plot_grid(  Lim1press,Lim2press, Lim3press, Lim1pulse,Lim2pulse,  Lim3pulse,Lim1pulsepress, Lim2pulsepress,Lim3pulsepress,  ncol = 3 )
corrRela

ggsave(plot = corrRela, file = here('output/simulations_corrRelatCon.png'), width = 12, height = 10)
## Absolute Contributions ##

# press 
Lim1RRpress <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1RRpress

Lim2RRpress <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2RRpress

Lim3RRpress <- stab.auc %>%
  filter(Model == 'press' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3RRpress


# pulse 
Lim1RRpulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1RRpulse

Lim2RRpulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2RRpulse

Lim3RRpulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3RRpulse


# pulsepress 
Lim1RRpulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1RRpulsepress

Lim2RRpulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit2') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim2RRpulsepress

Lim3RRpulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit3') %>%
  filter(species == 'species1' )%>%
  ggscatter(.,  x = 'mean.con.pi',y = 'AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Relative Contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim3RRpulsepress

corrAbs<- plot_grid(ncol = 3, Lim1press,Lim2press,Lim3press, Lim1pulse, Lim2pulse,Lim3pulse,  Lim1pulsepress,  Lim2pulsepress,Lim3pulsepress)
corrAbs

ggsave(plot = corrAbs, file = here('output/simulations_corrAbsCon.png'), width = 12, height = 10)


#### absolute Values for correlations ####

## Magnitude of relative contributions ##
# press 
Lim1press <- stab.auc %>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  filter(Model == 'press' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1press

Lim23press <- stab.auc %>%
  filter(Model == 'press' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23press


# pulse 
Lim1pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit1') %>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1pulse

Lim23pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23pulse

# pulsepress 
Lim1pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit1') %>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim1pulsepress

Lim23pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.pi = abs(AUC.pi)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.pi', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of relative contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23pulsepress


corrRelaPos<- plot_grid(  Lim1press,Lim23press, Lim1pulse,Lim23pulse,Lim1pulsepress, Lim23pulsepress,ncol = 2 )
corrRelaPos

ggsave(plot = corrRelaPos, file = here('output/simulations_corrRelatConAbsValues.png'), width = 8, height = 10)

## Magnitude of absolute contributions ##

# press 
Lim1press <- stab.auc %>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  filter(Model == 'press' & Limit == 'Limit1') %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09, label.y = 0.55)           # Add correlation coefficient
Lim1press

Lim23press <- stab.auc %>%
  filter(Model == 'press' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23press


# pulse 
Lim1pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit == 'Limit1') %>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09,  label.y = 0.65)           # Add correlation coefficient
Lim1pulse

Lim23pulse <- stab.auc %>%
  filter(Model == 'pulse' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23pulse

# pulsepress 
Lim1pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit == 'Limit1') %>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09, label.y = 1.3)           # Add correlation coefficient
Lim1pulsepress

Lim23pulsepress <- stab.auc %>%
  filter(Model == 'pulsepress' & Limit != 'Limit1') %>%
  filter(species == 'species1' )%>%
  mutate(pos.AUC.RR = abs(AUC.RR)) %>%
  ggscatter(.,  x = 'mean.con.pi',y = 'pos.AUC.RR', add = 'reg.line',alpha = 0.6, cor.method = 'spearmen',conf.int = T, ylab = 'Magnitude of absolute contribution', xlab = 'Relative dominance')  +
  stat_cor( label.x = 0.09)           # Add correlation coefficient
Lim23pulsepress


corrAbsPos<- plot_grid(  Lim1press,Lim23press, Lim1pulse,Lim23pulse,Lim1pulsepress, Lim23pulsepress,ncol = 2 )
corrAbsPos

ggsave(plot = corrAbsPos, file = here('output/simulations_corrAbsConAbsValues.png'), width = 8, height = 10)



### dominance ~competitiveness ###
ggscatter(stab.auc,  y = 'mean.con.pi',x = 'inv_relAlpha', add = 'reg.line', conf.int = T, ylab = 'Relat. competitiveness',xlab = 'Relat. dominance')  +
  stat_cor( label.x = 1)           # Add correlation coefficient
#ggsave(plot = last_plot(), file = here('output/Supplement_correlationDominanceAlpha.png'), width = 5, height = 4 )


#### Figure 3 dominance ~AUC  ####

# adjust decimal points:
scaleFUN <- function(x) sprintf("%.2f", x)

str(stab.auc)
rectLim1R <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                      x1 =c(0.044,NA, NA) ,
                      x2 = c(0.354,NA, NA) ,
                      y1 = c(0.85,NA, NA) ,
                      y2 = c(1.12,NA, NA) )


TextLim1rect  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                               x = 0.05,
                               y = 1,
                               text = c('All equally sensitive','',''))

Dom1R<-ggplot(subset(stab.auc, Limit == 'Limit1'), aes(x=mean.con.pi, y=AUC.RR,
                                                           col=species)) +
    geom_hline(yintercept=0, col="grey")+
    geom_point(alpha= .3, size = 2)+
    geom_text(data = TextLim1rect, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
    geom_rect(data = rectLim1R, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
    scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
    labs(y = 'Absolute contribution to stability', x = "")+
    facet_wrap(~Model, ncol = 3)+
    scale_y_continuous(limits = c(-1.7, 1.2), breaks = c(-1.5,-1,-0.5,0,0.5,1),labels=scaleFUN)+
    theme_bw() +
    theme(strip.background =element_rect(),
          strip.text = element_text(size = 12, face = 'bold'))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(text = element_text(size = 16)) +
    theme(legend.position="none")
Dom1R
  
rectLim1pi <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                    x1 =c(0.044,NA, NA) ,
                    x2 = c(0.354,NA, NA) ,
                    y1 = c(0.085,NA, NA) ,
                    y2 = c(0.105,NA, NA) )

TextLim1rectPi  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                        x = 0.05,
                        y = 0.095,
                        text = c('All equally sensitive','',''))

Dom1pi<-ggplot(subset(stab.auc, Limit == 'Limit1'), aes(x=mean.con.pi, y=AUC.pi,
                                                            col=species)) +
    geom_hline(yintercept=0, col="grey")+
    geom_point(alpha= .3, size = 2)+
  #  geom_text(data = TextLim1rectPi, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
  #  geom_rect(data = rectLim1pi, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
    scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
    scale_y_continuous(limits = c(-0.1, 0.11), breaks = c(-0.1,-0.05,0,0.05,0.1),labels=scaleFUN)+
    labs(y = 'Relative contribution to stability', x = " ")+
    facet_wrap(~Model, ncol = 3)+
    theme_bw() +
    theme(strip.background =element_rect(),
          strip.text = element_text(size = 12, face = 'bold'))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(text = element_text(size = 16)) +
    theme(legend.position="none")
Dom1pi
  

# Limit 2 #
rectLim2R <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                    x1 =c(0.044,NA, NA) ,
                    x2 = c(0.354,NA, NA) ,
                    y1 = c(26.8,NA, NA) ,
                    y2 = c(32.89,NA, NA) )


TextLim2rect  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                        x = 0.05,
                        y = 30,
                        text = c('One tolerant sp.','',''))

Dom2R<-ggplot(subset(stab.auc, Limit == 'Limit2'), aes(x=mean.con.pi, y=AUC.RR,
                                                           col=species)) +
    scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
    scale_y_continuous(limits = c(-34, 34), breaks = c(-30,-15,0,15,30),labels=scaleFUN)+
    geom_hline(yintercept=0, col="grey")+
    geom_point(alpha= .3, size = 2)+
    geom_text(data = TextLim2rect, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
    geom_rect(data = rectLim2R, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
    labs(y = 'Absolute contribution to stability', x = " ")+
    facet_wrap(~Model, ncol = 3)+
    theme_bw() +
    theme(strip.text.x = element_text(color = 'white'),
          strip.background =element_rect(fill = 'white', linetype = 0))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(text = element_text(size = 16)) +
    theme(legend.position="none")

Dom2R
  

rectLim2pi <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                    x1 =c(0.044,NA, NA) ,
                    x2 = c(0.354,NA, NA) ,
                    y1 = c(3.6,NA, NA) ,
                    y2 = c(4.35,NA, NA) )


TextLim2rectPi  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                        x = 0.05,
                        y = 4,
                        text = c('One tolerant sp.','',''))

Dom2pi<-ggplot(subset(stab.auc, Limit == 'Limit2'), aes(x=mean.con.pi, y=AUC.pi,
                                                            col=species)) +
    scale_color_manual(values = c('#68789E', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
    scale_y_continuous(limits = c(-4.5, 4.5), breaks = c(-4,-2,0,2,4),labels=scaleFUN)+
    geom_hline(yintercept=0, col="grey")+
    geom_point(alpha= .3, size = 2)+
 # geom_text(data = TextLim2rectPi, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
 # geom_rect(data = rectLim2pi, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
  labs(y = 'Relative contribution to stability', x = "")+
    facet_wrap(~Model, ncol = 3)+
    theme_bw() +
    theme(strip.text.x = element_text(color = 'white'),
          strip.background =element_rect(fill = 'white', linetype = 0))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(text = element_text(size = 16)) +
    theme(legend.position="none")
Dom2pi
  

# Limit 3 #
rectLim3R <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                    x1 =c(0.044,NA, NA) ,
                    x2 = c(0.354,NA, NA) ,
                    y1 = c(17.10,NA, NA) ,
                    y2 = c(22.69,NA, NA) )


TextLim3rect  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                        x = 0.05,
                        y = 20,
                        text = c('One sensitive sp.','',''))

Dom3R<-ggplot(subset(stab.auc, Limit == 'Limit3'), aes(x=mean.con.pi, y=AUC.RR,
                                                           col=species)) +
    geom_hline(yintercept=0, col="grey")+
    geom_point(alpha= .3, size = 2)+
    geom_text(data = TextLim3rect, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
    geom_rect(data = rectLim3R, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
    scale_color_manual(values = c('#D56060', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
    scale_y_continuous(limits = c(-30, 25), breaks = c(-30,-20,-10,0,10,20),labels=scaleFUN)+
    labs(y='Absolute contribution to stability', x = "Relative dominance")+
    facet_wrap(~Model, ncol = 3)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(strip.text.x = element_text(color = 'white'),
          strip.background =element_rect(fill = 'white', linetype = 0))+
    theme(text = element_text(size = 16)) +  theme(legend.position="none")
  
Dom3R

rectLim3pi <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                    x1 =c(0.044,NA, NA) ,
                    x2 = c(0.354,NA, NA) ,
                    y1 = c(1.710,NA, NA) ,
                    y2 = c(2.269,NA, NA) )


TextLim3rectPi  <- tibble(Model = c('press', 'pulse', 'pulsepress'),
                        x = 0.05,
                        y = 2,
                        text = c('One sensitive sp.','',''))

Dom3pi<-ggplot(subset(stab.auc, Limit == 'Limit3'), aes(x=mean.con.pi, y=AUC.pi,
                                                            col=species)) +
  geom_hline(yintercept=0, col="grey")+
  geom_point(alpha= .3, size = 2)+
 # geom_text(data = TextLim3rectPi, inherit.aes = F, aes(x = x, y = y, label = text), size = 4, hjust = 0) +
 # geom_rect(data = rectLim3pi, inherit.aes = F,mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, colour = '#bababa') +
  scale_color_manual(values = c('#D56060', '#BEBEBE', '#BEBEBE', '#BEBEBE', '#BEBEBE'))+
  scale_y_continuous(limits = c(-3, 2.5), breaks = c(-3.0,-2.0,-1.0,0,1.0,2.0),labels=scaleFUN)+
    labs(y= 'Relative contribution to stability', x = "Relative dominance")+
    facet_wrap(~Model, ncol = 3)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
    theme(strip.text.x = element_text(color = 'white'),
          strip.background =element_rect(fill = 'white', linetype = 0))+
    theme(text = element_text(size = 16)) +  theme(legend.position="")
  
Dom3pi#ggsave(plot = RRpi, 'RRPi.png',width = 8, height = 4)

par(mar=c (5.1, 4.1, 4.1, 2.1))
plot_grid(Dom1R,Dom1pi, Dom2R,Dom2pi,Dom3R,Dom3pi,hjust = -1, ncol = 2, nrow = 3, labels = c( '(a)',' ', '(b)',' ', '(c)'), rel_widths = c(1,1))
ggsave(plot = last_plot(), width = 12, height = 12, file = here('output/Fig.3_Dominance_allDist.png'))
