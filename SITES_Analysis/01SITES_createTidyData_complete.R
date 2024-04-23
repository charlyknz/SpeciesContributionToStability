#### R script to calculate species contributions to stability for SITES Data ####

#packages
library(tidyverse)
library(here)
library(MESS)

#### import data ####
zoo <- read.csv2('SITES_Data/zooplankton.csv', sep = ';')
str(zoo)
names(zoo)

#check data structure 
dum<-zoo %>% 
  filter(Treatment != 'Lake' & Replicate != 'lake') %>%# reduce dataset to experimental dataset
  group_by(Lake, Experiment, Treatment, Exp_day, Taxa) %>% 
  summarise(replicates = n())
print(dum)
rm(dum)#remove 

### Adjust data ###
# to make sure our datasets from all different stations have reported species in the same way, 
# we use complete grouped by lake to fill up with missing species - we do this by species to make sure we only include species which could geographically occur.
# these NAs are then set to 0 to calculate the correct mean value 
# example: In 1 of 4 control replicates we find 4 daphnia: mean(Na,Na,Na,4) = 4 whereas mean(0,0,0,4) = 1

#use cleaned biomass values & add rows for species which had not been included because of Biomass = = 
zoo.data <- zoo %>%
  filter(Treatment != 'Lake' ) %>% # reduce dataset to experimental dataset
  filter(!Replicate %in%c('Lake', 'lake')) %>%
  select(Lake, Experiment, Treatment, Mean_lenght, Exp_day, Replicate, Taxa, Abundance,Clean_Biomass)%>%
  group_by(Lake, Experiment)%>%
  complete(Treatment, Exp_day,Taxa, Replicate, fill = list(Biomass = 0)) %>% # have Biomass which is NA as 0 to calculate the correct mean here 
  mutate(Clean_Biomass = as.numeric(Clean_Biomass)) %>%
  rename(Biomass = Clean_Biomass)
zoo.data$Biomass[is.na(zoo.data$Biomass)]<-0
summary(zoo.data)

### check if complete() worked ###
which(is.na(zoo.data$Biomass))

#saveRDS(zoo.data, file = here('SITES_Data/ZooData.RDS'))

#### Calculation of Response ratios ####

## Control data ##
# calculate a mean response for control to calculate RR
control.zoo <- zoo.data %>%
  filter(Treatment == 'C') %>%
  group_by(Lake, Experiment, Taxa, Exp_day) %>%
  summarise(con.bio = mean(Biomass, na.rm = T))


## merge control and treatment data & calculate overall response ##
zoo.data2 <- zoo.data %>%
  filter(Treatment != 'C') %>% #remove control data from Treatment data 
  merge(., control.zoo, by = c('Lake', 'Experiment', 'Exp_day', 'Taxa')) %>%
  group_by(Lake, Experiment, Treatment, Exp_day, Replicate) %>%
  mutate(treat.tot = sum(Biomass, na.rm = T),
         con.tot = sum(con.bio, na.rm = T),
         deltabm.tot = (treat.tot - con.tot)/(treat.tot+con.tot))  #calculate tot.RR
hist(zoo.data2$deltabm.tot)


# remove rows where treatment and control are both 0
zoo.data2$dummyRR <- zoo.data2$con.bio + zoo.data2$Biomass
zoo.data2 <- dplyr::filter(zoo.data2, dummyRR!=0)

## calculate species-specific measures ##
zoo.data3 <- zoo.data2 %>%
  select(-dummyRR)  %>%
  mutate(USI = paste(Lake, Experiment, Taxa, Replicate, Treatment, sep = '_')) #create Unique Identifier


# create species specific change in biomass  (RR) 
zoo.data3$RR<-(zoo.data3$Biomass-zoo.data3$con.bio)/(zoo.data3$Biomass+zoo.data3$con.bio)

# create species specific contribution (pi) to biomass for control and treat  (treat / tot)         ## contribution 
zoo.data3$treat.pi <- zoo.data3$Biomass/zoo.data3$treat.tot
zoo.data3$con.pi <- zoo.data3$con.bio/zoo.data3$con.tot

# create species specific LRR and difference for pi           (delta pi)            
zoo.data3$delta.pi <- zoo.data3$treat.pi-zoo.data3$con.pi


#### AUC calculation ####
names(zoo.data3)

zoo.data3 <- zoo.data3[order(zoo.data3$Exp_day),]

#create an empty data frame (then our loop is faster)
zoo.stab.auc <- data.frame()


# the following loops cycle through all unique cases (USI)
USI <- unique(zoo.data3$USI)


for(i in 1:length(USI)){
  temp<-zoo.data3[zoo.data3$USI==USI[i], ]#creates a temporary data frame for each case
  if(dim(temp)[1]>2){#does the next step only if at least 3 data points are present
    AUC.RR<-auc(temp$Exp_day, temp$RR,  from = min(temp$Exp_day, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                type = c("linear"),absolutearea = FALSE)
    AUC.pi<-auc(temp$Exp_day, temp$delta.pi, from = min(temp$Exp_day, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                type = c("linear"),absolutearea = FALSE)
    AUC.totRR<-auc(temp$Exp_day, temp$deltabm.tot, from = min(temp$Exp_day, na.rm = T), to = max(temp$Exp_day, na.rm = TRUE),
                   type = c("linear"),absolutearea = FALSE)
    mean.con.pi <- mean(temp$con.pi, na.rm = T)
    mean.treat.pi <- mean(temp$treat.pi, na.rm = T)
    zoo.stab.auc<-rbind(zoo.stab.auc,
                        data.frame(temp[1,c(1,2,4,5,6,13)],
                                   AUC.RR,
                                   AUC.pi,
                                   AUC.totRR,
                                   mean.con.pi,
                                   mean.treat.pi))
    rm(temp)
  }
}

str(zoo.stab.auc)


# remove undefined Taxa, pelagic species and Mosquito larvae #
unique(zoo.stab.auc$Taxa)

zoo.stab.auc.prefin <- zoo.stab.auc %>% 
  filter(!Taxa %in% c("Nauplii", "Small cladoceran", 'Sida', 'Chaoborus'))

## adjust Lake and Disturbance Type label ##
unique(zoo.stab.auc.prefin$Treatment)
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

### write csv ###
write_csv(zoo.stab.auc.prefin, file = here('OutputSubmission/AUCdata_3.csv'))


#### AUC Plots ####
ggplot(zoo.stab.auc.prefin,  aes(AUC.pi, AUC.RR,color = Taxa )) +
  geom_point(aes(shape = as.factor(Experiment)), alpha = 0.5) +
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

## Fig. 4 - Mean absolute and relative contributions ##
# calculate mean +- SE 
zoo.stab.auc.mean <- zoo.stab.auc.prefin %>% 
  group_by(Taxa) %>%
  mutate(mean.dom = mean(mean.con.pi)) %>% #calculate relat. mean dominance for point size
  dplyr::group_by(Taxa, Lake,Experiment, Treatment, mean.dom) %>% 
  dplyr::summarise(mean.AUC.RR = mean(AUC.RR, na.rm = T),
                mean.AUC.pi = mean(AUC.pi, na.rm = T),
                sd.AUC.RR = sd(AUC.RR, na.rm = T),
                sd.AUC.pi = sd(AUC.pi, na.rm = T),
                se.AUC.RR = sd.AUC.RR/sqrt(n()),
                se.AUC.pi = sd.AUC.pi/sqrt(n()))  

unique(zoo.stab.auc.mean$Taxa)

# plot 
ggplot(zoo.stab.auc.mean,aes(mean.AUC.pi, mean.AUC.RR,color = Taxa )) +
  geom_errorbarh(aes(xmax = mean.AUC.pi + se.AUC.pi, xmin = mean.AUC.pi - se.AUC.pi))+
  geom_errorbar(aes(ymax = mean.AUC.RR + se.AUC.RR, ymin = mean.AUC.RR - se.AUC.RR))+
  geom_point(aes(shape = as.factor(Experiment)), alpha = 0.5, size = 2.5) +
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
ggsave(plot = last_plot(), file = here('OutputSubmission/Fig4.tiff'),width = 14, height = 10)

### look at two species exemplary ####
zoo.stab.auc.mean %>%
  filter(Taxa %in%c('Bosmina', 'Keratella') ) %>%
  ggplot(.,aes(mean.AUC.pi, mean.AUC.RR,color = Taxa )) +
  geom_point(aes(shape = as.factor(Experiment)), alpha = 0.5, size = 2) +
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


#### Sector counts ####
### calculate how many percent of taxa are present in each sector 

zoo.stab.auc.mean$Sector<- NA
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi>0&zoo.stab.auc.mean$mean.AUC.RR>0]<-1
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi>0&zoo.stab.auc.mean$mean.AUC.RR<0]<-2
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi<0&zoo.stab.auc.mean$mean.AUC.RR<0]<-3
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi<0&zoo.stab.auc.mean$mean.AUC.RR>0]<-4
zoo.stab.auc.mean$Sector[zoo.stab.auc.mean$mean.AUC.pi== 0 & zoo.stab.auc.mean$mean.AUC.RR==0] <-0

# create df 
sector.count<- zoo.stab.auc.mean%>%
  group_by( Sector) %>%
  count()%>%
  rename(N = n) %>%
  drop_na(Sector) %>%
  ungroup()%>%
  mutate(sum = sum(N, na.rm = T),
         relN = (N/sum)*100)

sector.count # output

#### Compensatory Dynamics analysis ####
names(zoo.stab.auc.prefin)


zoo.stab.auc.prefin$Sector<- NA
zoo.stab.auc.prefin$Sector[zoo.stab.auc.prefin$AUC.pi>0&zoo.stab.auc.prefin$AUC.RR>0]<-1
zoo.stab.auc.prefin$Sector[zoo.stab.auc.prefin$AUC.pi>0&zoo.stab.auc.prefin$AUC.RR<0]<-2
zoo.stab.auc.prefin$Sector[zoo.stab.auc.prefin$AUC.pi<0&zoo.stab.auc.prefin$AUC.RR<0]<-3
zoo.stab.auc.prefin$Sector[zoo.stab.auc.prefin$AUC.pi<0&zoo.stab.auc.prefin$AUC.RR>0]<-4
zoo.stab.auc.prefin$Sector[zoo.stab.auc.prefin$AUC.pi== 0 & zoo.stab.auc.prefin$AUC.RR==0] <-0

AUCsum <- zoo.stab.auc.prefin %>%
  group_by(Lake, Experiment, Treatment, Sector, Replicate) %>%
  summarise(sum.RR = sum(AUC.RR),
            sum.pi = sum(AUC.pi))


spring <- AUCsum %>%
  filter(Experiment == 'spring') %>%
  group_by(Lake, Experiment, Treatment, Sector) %>%
  summarise(mean.RR = mean(sum.RR),
            mean.pi = mean(sum.pi),
            sd.RR = sd(sum.RR),
            sd.pi = sd(sum.pi),
            se.RR = sd.RR / sqrt(n()),
            se.pi = sd.pi / sqrt(n())) %>%
ggplot(., aes(x = Sector, y = mean.pi, fill = as.factor(Sector)))+
  geom_errorbar(aes(ymin = mean.pi - se.pi, ymax = mean.pi + se.pi))+
  geom_col(color = 'black')+
  scale_fill_manual(values = c('#68789E', '#BEBEBE','#D56060', 'white'))+
  labs(y = 'Total Relative Contribution to Stability', fill = 'Sector', title = 'Spring Experiment')+
  geom_hline(yintercept = 0)+
  facet_grid(~Lake~Treatment)+
  theme_bw()

spring.rr <- AUCsum %>%
  summarise(mean.RR = mean(sum.RR),
            mean.pi = mean(sum.pi),
            sd.RR = sd(sum.RR),
            sd.pi = sd(sum.pi),
            se.RR = sd.RR / sqrt(n()),
            se.pi = sd.pi / sqrt(n())) %>%
  filter(Experiment == 'spring') %>%
  ggplot(., aes(x = Sector, y = mean.RR, fill = as.factor(Sector)))+
  geom_errorbar(aes(ymin = mean.RR - se.RR, ymax = mean.RR + se.RR))+
  geom_col(color = 'black')+
  scale_fill_manual(values = c('#68789E', '#BEBEBE','#D56060', 'white'))+
  labs(y = 'Total Absolute Contribution to Stability', fill = 'Sector', title = 'Spring Experiment')+
  geom_hline(yintercept = 0)+
  facet_grid(~Lake~Treatment)+
  theme_bw()


cowplot::plot_grid(spring, spring.rr, ncol = 1, labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = here('OutputSubmission/FigS3_AUCcontribution_spring.tiff'), width = 7, height = 14)

summer <- AUCsum %>%
  filter(Experiment == 'summer') %>%
  group_by(Lake, Experiment, Treatment, Sector) %>%
  summarise(mean.RR = mean(sum.RR),
            mean.pi = mean(sum.pi),
            sd.RR = sd(sum.RR),
            sd.pi = sd(sum.pi),
            se.RR = sd.RR / sqrt(n()),
            se.pi = sd.pi / sqrt(n())) %>%
  drop_na(Sector)%>%
  ggplot(., aes(x = Sector, y = mean.pi, fill = as.factor(Sector)))+
  geom_errorbar(aes(ymin = mean.pi - se.pi, ymax = mean.pi + se.pi))+
  geom_col(color = 'black')+
  scale_fill_manual(values = c('#68789E', '#BEBEBE','#D56060', 'white'))+
  labs(y = 'Total Relative Contribution to Stability', fill = 'Sector', title = 'Summer Experiment')+
  geom_hline(yintercept = 0)+
  facet_grid(~Lake~Treatment)+
  theme_bw()

summer.rr <- AUCsum %>%
  summarise(mean.RR = mean(sum.RR),
            mean.pi = mean(sum.pi),
            sd.RR = sd(sum.RR),
            sd.pi = sd(sum.pi),
            se.RR = sd.RR / sqrt(n()),
            se.pi = sd.pi / sqrt(n())) %>%
  filter(Experiment == 'summer') %>%
  drop_na(Sector)%>%
  ggplot(., aes(x = Sector, y = mean.RR, fill = as.factor(Sector)))+
  geom_errorbar(aes(ymin = mean.RR - se.RR, ymax = mean.RR + se.RR))+
  geom_col(color = 'black')+
  scale_fill_manual(values = c('#68789E', '#BEBEBE','#D56060', 'white'))+
  labs(y = 'Total Absolute Contribution to Stability', fill = 'Sector', title = 'Summer Experiment')+
  geom_hline(yintercept = 0)+
  facet_grid(~Lake~Treatment)+
  theme_bw()


cowplot::plot_grid(summer, summer.rr, ncol = 1, labels = c('(a)', '(b)'))
ggsave(plot = last_plot(), file = here('OutputSubmission/FigS4_AUCcontribution_summer.tiff'), width = 7, height = 14)


### Total biomass over time ###
Total.zoo <- zoo.data 

Total.zoo[Total.zoo$Lake == "Asa", "Lake"] <- "Feresjörn"
Total.zoo[Total.zoo$Lake == "Erk", "Lake"] <- "Erken"
Total.zoo[Total.zoo$Lake == "Sva", "Lake"] <- "Stortjärn"
Total.zoo[Total.zoo$Lake == "Sko", "Lake"] <- "Erssjön"
Total.zoo[Total.zoo$Lake == "Bol", "Lake"] <- "Bolmen"

Total.zoo[Total.zoo$Treatment == "C", "Treatment"] <- "Control"
Total.zoo[Total.zoo$Treatment == "F", "Treatment"] <- "Pulse"
Total.zoo[Total.zoo$Treatment == "FS", "Treatment"] <- "Pulse & Press"
Total.zoo[Total.zoo$Treatment == "S", "Treatment"] <- "Press"

Total.zoo$Experiment[Total.zoo$Experiment == 1] <- 'spring'
Total.zoo$Experiment[Total.zoo$Experiment == 2] <- 'summer'

Total.zoo%>%
  group_by(Lake, Experiment, Treatment, Exp_day, Replicate)%>%
  summarise(sum = sum(Biomass))%>%
  group_by(Lake, Experiment, Treatment, Exp_day)%>%
  summarise(mean = mean(sum),
            sd = sd(sum),
            se = sd/sqrt(n())) %>%
  ggplot(., aes (x = Exp_day, y = mean, color = Treatment))+
  geom_line()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .7)+
  geom_point()+
  scale_color_manual(values = c('black', '#DF536B', '#61D04F', '#2297E6'))+
  labs(x= 'Time [days]', y = expression(Total~Biomass~'['~mu*g*L^{-1}~']'))+
  facet_grid(~Lake~Experiment, scale = 'free_y')+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  

ggsave(plot = last_plot(), file = here('OutputSubmission/FigS5_TotalBiomass.tiff'), width = 8, height = 8)


