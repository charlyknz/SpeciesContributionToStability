#### R script to calculate AUCs for SITES Data ####
# 19.05.2023

library(tidyverse)
library(here)
library(MESS)

#### import data ####
zoo <- read.csv2('SITES_Data/zooplankton.csv', sep = ';')
str(zoo)

#remove 

#use cleaned biomass values
zoo.data <- zoo %>%
  select(-Biomass) %>%
  mutate(Clean_Biomass = as.numeric(Clean_Biomass)) %>%
  rename(Biomass = Clean_Biomass) %>%
  filter(Treatment != 'Lake') # reduce dataset to experimental dataset

zoo.data$Biomass[zoo.data$Biomass == 0 & is.na(zoo.data$Abundance)] <- NA

levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Asa"] <- "Fer"    # Feresjoen
levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Sko"] <- "Ers"    # Erssjoen
levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Sva"] <- "Sto"    # Storjaen


## Control data ##
control.zoo <- zoo.data %>%
  filter(Treatment == 'C') %>%
  group_by(Lake, Experiment, Taxa, Sampling) %>%
  summarise(con.bio = mean(Biomass, na.rm = T))


## merge control and treatment data & calculate overall response ##
zoo.data2 <- zoo.data %>%
  filter(Treatment != 'C') %>% #remove control data from Treatment data 
  merge(., control.zoo, by = c('Lake', 'Experiment', 'Sampling', 'Taxa')) %>%
  group_by(Lake, Experiment, Treatment, Sampling, Replicate) %>%
  mutate(treat.tot = sum(Biomass, na.rm = T),
         con.tot = sum(con.bio, na.rm = T),
         LRR.tot = log(treat.tot/con.tot),
         deltabm.tot = (treat.tot - con.tot)/(treat.tot+con.tot))  #calculate tot.RR
hist(zoo.data2$deltabm.tot)


# remove rows where treatment and control are both 0
zoo.data2$dummyRR <- zoo.data2$con.bio + zoo.data2$Biomass
zoo.data2 <- dplyr::filter(zoo.data2, dummyRR!=0)

## calculate species-specific measures ##
zoo.data3 <- zoo.data2 %>%
  select(-dummyRR)  %>%
  mutate(USI = paste(Lake, Experiment, Taxa, Replicate, Treatment, sep = '_')) #create Unique Identifier


# create species specific Log Response Ratio (LRR) for biomass
zoo.data3$LRR<-log(zoo.data3$Biomass/zoo.data3$con.bio)


# create species specific change in biomass  (RR) 
zoo.data3$RR<-(zoo.data3$Biomass-zoo.data3$con.bio)/(zoo.data3$Biomass+zoo.data3$con.bio)

# the absence of species in control or treatment creates INF, get rid of these
zoo.data3$LRR[zoo.data3$LRR=="Inf"]<-NA
zoo.data3$LRR[zoo.data3$LRR=="-Inf"]<-NA

# create species specific contribution (pi) to biomass for control and treat  (treat / tot)         ## contribution 
zoo.data3$treat.pi <- zoo.data3$Biomass/zoo.data3$treat.tot
zoo.data3$con.pi <- zoo.data3$con.bio/zoo.data3$con.tot

# create species specific LRR and difference for pi           (delta pi)            
zoo.data3$delta.pi <- zoo.data3$treat.pi-zoo.data3$con.pi


# create the deviance between species and community effect sizes                     
zoo.data3$LRR.diff <- zoo.data3$LRR-zoo.data3$LRR.tot                                


#### AUC calculation ####
names(zoo.data3)

#order data by day to start with 0 
zoo.data3 <- zoo.data3[order(zoo.data3$Exp_day),]

#create an empty data frame (then our loop is faster)
zoo.stab.auc <- data.frame()


# the following loops cycle through all unique cases (USI)
USI <- unique(zoo.data3$USI)

for(i in 1:length(USI)){
  temp<-zoo.data3[zoo.data3$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$Exp_day, temp$RR,  from = min(0, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.pi<-auc(temp$Exp_day, temp$delta.pi, from = min(0, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                type = c("spline"),absolutearea = FALSE)
    AUC.totRR<-auc(temp$Exp_day, temp$deltabm.tot, from = min(0, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                   type = c("spline"),absolutearea = FALSE)
    mean.con.pi <- mean(temp$con.pi, na.rm = T)
    zoo.stab.auc<-rbind(zoo.stab.auc,
                        data.frame(temp[1,c(1,2,4,9,10,11)],
                                   AUC.RR,
                                   AUC.pi,
                                   AUC.totRR,
                                   mean.con.pi))
    rm(temp)
  }
}

str(zoo.stab.auc)


# remove undefined Taxa #
unique(zoo.stab.auc$Taxa)

zoo.stab.auc.prefin <- zoo.stab.auc %>% 
  filter(!Taxa %in% c("Nauplii", "Small cladoceran"))

## adjust Lake and Disturbance Type label ##
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Asa", "Lake"] <- "Feresjörn"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Erk", "Lake"] <- "Erken"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Sva", "Lake"] <- "Stortjärn"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Sko", "Lake"] <- "Erssjön"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Bol", "Lake"] <- "Bolmen"

zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "F", "Treatment"] <- "Pulse"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "FS", "Treatment"] <- "Pulse & Press"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "S", "Treatment"] <- "Press"

zoo.stab.auc.prefin$Experiment[zoo.stab.auc.prefin$Experiment == 1] <- 'spring'
zoo.stab.auc.prefin$Experiment[zoo.stab.auc.prefin$Experiment == 2] <- 'summer'


#write.csv2(zoo.stab.auc.prefin, 'AUCdata_2.csv')


#### AUC Plots ####
ggplot(zoo.stab.auc.prefin,  aes(AUC.pi, AUC.RR,color = Taxa )) +
  geom_point(aes(shape = as.factor(Experiment),size = mean.con.pi), alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.25) +                                      
  geom_hline(yintercept = 0, alpha = 0.25) +
  facet_grid(Lake~Treatment, scales = 'free') +
  labs(x = expression(AUC.~Delta~'pi'~'('~Compositional~Stability~')'), y = "Functional Stability (AUC.RR)", shape = 'Exp') +  
  theme_bw()+
  theme(legend.position = 'right')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(plot.title = element_text(size=14),
        plot.subtitle = element_text(size=13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))
#ggsave(plot = last_plot(), file = here('output/ZooAUC_reps.png'), width = 10, height = 8)

## Fig. 4 - AUC.RR and pi with errorbars ##
zoo.stab.auc.mean <- zoo.stab.auc.prefin %>% 
  group_by(Taxa) %>%
  mutate(mean.dom = mean(mean.con.pi)) %>% #calculate mean dominance for point size
  dplyr::group_by(Taxa, Lake,Experiment, Treatment, mean.dom) %>% # for single dots use mutate here!! 
  dplyr::summarise(mean.AUC.RR = mean(AUC.RR, na.rm = T),
                mean.AUC.pi = mean(AUC.pi, na.rm = T),
                sd.AUC.RR = sd(AUC.RR, na.rm = T),
                sd.AUC.pi = sd(AUC.pi, na.rm = T),
                se.AUC.RR = sd.AUC.RR/sqrt(n()),
                se.AUC.pi = sd.AUC.pi/sqrt(n()))  

ggplot(zoo.stab.auc.mean,aes(mean.AUC.pi, mean.AUC.RR,color = Taxa )) +
  geom_point(aes(shape = as.factor(Experiment), size = mean.dom), alpha = 0.5) +
  geom_errorbarh(aes(xmax = mean.AUC.pi + se.AUC.pi, xmin = mean.AUC.pi - se.AUC.pi))+
  geom_errorbar(aes(ymax = mean.AUC.RR + se.AUC.RR, ymin = mean.AUC.RR - se.AUC.RR))+
  geom_vline(xintercept = 0, alpha = 0.25) +                                      
  geom_hline(yintercept = 0, alpha = 0.25) +
  facet_grid(Lake~Treatment, scales = 'free') +
  labs(y = 'Absolute contribution to stability', x = "Relative contribution to stability", size = 'Dominance', shape = 'Season') +  
  theme_bw()+
  theme(legend.position = 'right')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 17,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12, face = "plain", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 17, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12, face = "plain", colour = "black", angle = 0, hjust = 0.4)) +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 15, face = 'bold'),
        strip.text.y  = element_text(size = 15, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
ggsave(plot = last_plot(), file = here('output/Fig4.png'),width = 12, height = 8)


#### Sector counts ####
zoo.stab.auc.mean$Sector<- NA
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi>0&zoo.stab.auc.mean$mean.AUC.RR>0]<-1
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi>0&zoo.stab.auc.mean$mean.AUC.RR<0]<-2
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi<0&zoo.stab.auc.mean$mean.AUC.RR<0]<-3
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi<0&zoo.stab.auc.mean$mean.AUC.RR>0]<-4
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi== 0 & zoo.stab.auc.mean$mean.AUC.RR==0] <-0


sector.count<- zoo.stab.auc.mean%>%
  group_by( Sector) %>%
  count()%>%
  rename(N = n) %>%
  drop_na(Sector) %>%
  ungroup()%>%
  mutate(sum = sum(N, na.rm = T),
         relN = (N/sum)*100)
