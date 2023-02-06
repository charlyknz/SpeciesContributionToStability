### Code for analyzing species stability
rm(list=ls())
graphics.off()

# the following libraries are needed
# please add all of them if they are not installed yet
library(cowplot)
library(tidyverse)
library(GGally)
library(MESS)
library(here)

#### load data ####

zoo <- read.csv2("SITES_Data/zooplankton.csv", sep=";") 
str(zoo)
names(zoo)



#### size data ####
zoo.size <- zoo %>%
  mutate(size = Mean_lenght*1000) %>% #convert size into um
  group_by(Zooplankton_group, Treatment, Taxa ) %>%
  summarise(
    mean.size = mean(size, na.rm = T),
    sd.size = sd(size, na.rm = T)) %>%
  arrange(Treatment, desc(mean.size)) %>%
  filter(!Taxa %in% c('small Cladoceran', 'Nauplii'))
write.csv(zoo.size, file = here('SITES_Data/zooplanktonMeanLength.csv'))

#### Biomass based Analysis ####
zoo <- zoo %>%
  dplyr::select( -Biomass) %>%
  dplyr::mutate(Clean_Biomass = as.numeric(Clean_Biomass)) %>%
  dplyr::rename(Biomass = Clean_Biomass)
zoo$Biomass[zoo$Biomass == 0 & is.na(zoo$Abundance)] <- NA

#some changes to the data
#reduce to experimental subset
zoo.data <- dplyr::filter(zoo, Treatment != 'Lake')

#make sure lakes are named correctly, always using first three letters for abbreviation

levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Asa"] <- "Fer"    # Feresjoen
levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Sko"] <- "Ers"    # Erssjoen
levels(zoo.data$Lake)[levels(zoo.data$Lake)=="Sva"] <- "Sto"    # Storjaen

#extract the control data and merge as an extra variable

con.zoo <- dplyr::filter(zoo.data, Treatment=="C")
con.zoo <-  dplyr::rename(con.zoo, "con.bio"= "Biomass")                                    # con.bio  -->  Biomass

#calculate mean biomass in control dataset 

con.zoo <- con.zoo %>%
  dplyr::group_by(Lake, Experiment, Sampling, Taxa) %>%
  dplyr::summarise(con.bio = mean(con.bio, na.rm = T))

#bringing in the control data to the main data

zoo.data2 <- merge(zoo.data,con.zoo, by = c("Lake","Experiment",                  # keeping data from the second (right) table and joins anything 
                                                        "Sampling","Taxa"))               # that matches from the first (left) table.
# take out the control rows

zoo.data2 <- dplyr::filter(zoo.data2, Treatment!="C")

# for comparison, calculate the overall response                                    ## data.tot -> Total Response on Forcing

zoo.data.tot<-zoo.data2%>% 
  dplyr::group_by(Lake, Experiment, Sampling, Treatment, Replicate) %>%
  dplyr:: summarise(treat.tot = sum(Biomass, na.rm = T),
                    con.tot = sum(con.bio, na.rm = T))

# LRR (Log Response Ratio) of community 

zoo.data.tot$LRR.tot<-log(zoo.data.tot$treat.tot/zoo.data.tot$con.tot)
zoo.data.tot$deltabm.tot<-(zoo.data.tot$treat.tot-zoo.data.tot$con.tot)/(zoo.data.tot$treat.tot+zoo.data.tot$con.tot)
hist(zoo.data.tot$LRR.tot)
hist(zoo.data.tot$deltabm.tot)

#create data set for species specific effect sizes
#take out those rows where biomass is 0 in both treatment (Biomass) + control (con.bio)

zoo.data3 <- zoo.data2 %>%
  dplyr::left_join(.,zoo.data.tot, by = c("Lake","Experiment","Sampling","Replicate","Treatment"))
zoo.data3$RR <- zoo.data3$Biomass+zoo.data3$con.bio
zoo.data3 <- dplyr::filter(zoo.data3, RR!=0)

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


# weight LRR  -  mean pi                                                            
zoo.data3$mean.pi <- 0.5*(zoo.data3$treat.pi+zoo.data3$con.pi)
zoo.data3$LRR.w <- zoo.data3$LRR*zoo.data3$mean.pi
zoo.data3$RR.w <- zoo.data3$RR*zoo.data3$mean.pi

# create the deviance between species and community effect sizes                     
zoo.data3$LRR.diff <- zoo.data3$LRR-zoo.data3$LRR.tot                                
zoo.data3$LRR.diff.w <- zoo.data3$LRR.diff*zoo.data3$mean.pi


# differentiate the most dominant species 

zoo.mean.dom <- zoo.data3 %>%
  dplyr:: group_by(Taxa) %>%
  dplyr:: summarise(mean=mean(mean.pi, na.rm = T))
zoo.mean.dom$dom <- zoo.mean.dom$mean  #*zoo.mean.dom$length                       
zoo.data3 <- merge(zoo.data3, dplyr::select(zoo.mean.dom, -mean), by = c("Taxa"))
hist(zoo.data3$dom)
str(zoo.data3)
# create mean and standard deviation for community LRR
zoo.tot.mn <- zoo.data.tot %>%
  dplyr::group_by(Lake, Experiment, Sampling, Treatment) %>%
  dplyr:: summarise(LRR.mean=mean(LRR.tot, na.rm = T),
                    LRR.sd=sd(LRR.tot, na.rm = T))
zoo.data3 <- merge(zoo.data3,zoo.tot.mn, by = c("Lake","Experiment","Sampling","Treatment"))
summary(zoo.data3)

unique(zoo.data3$Taxa)



####  AUC Metrics  ###
# collapse variation between replicates

zoo.data3.mn <- zoo.data3 %>%
  dplyr::group_by(Lake, Treatment, Experiment, Taxa, Exp_day, Sampling) %>%
  dplyr::summarise(across(c(Biomass:LRR.sd), mean))


#### AUC Loop ####

names(zoo.data3)
zoo.cols <- c("Lake","Experiment","Treatment","Taxa","Replicate")
zoo.data3 <- zoo.data3 %>% 
  dplyr::mutate(USI = paste(Lake, Experiment, Taxa, Replicate, Treatment, sep = "_"))     #create unique ID

#order data by day to start with 0 
zoo.data3 <- zoo.data3[order(zoo.data3$Exp_day),]
names(zoo.data3)
#create an empty data frame (then our loop is faster)
zoo.stab.auc <- data.frame()


# the following loops cycle through all unique cases (USI)

# Zooplankton
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
    AUC.LRR<-auc(temp$Exp_day, temp$LRR.mean, from = min(0, na.rm = TRUE), to = max(temp$Exp_day, na.rm = TRUE),
                   type = c("spline"),absolutearea = FALSE)
    mean.delta.pi<-mean(temp$delta.pi)
    mean.RR<-mean(temp$RR)
    mean.pi<-mean(temp$mean.pi)
    sd.delta.pi<-sd(temp$delta.pi)
    sd.RR<-sd(temp$RR)
    zoo.stab.auc<-rbind(zoo.stab.auc,
                        data.frame(temp[1,c(1,2,4,5,8,9,10,11,25,27,28,34)],
                                   AUC.RR,
                                   AUC.pi,
                                   AUC.totRR,
                                   AUC.LRR,
                                   mean.delta.pi,
                                   mean.RR,
                                   mean.pi,
                                   sd.delta.pi,
                                   sd.RR))
    rm(temp)
  }
}

summary(zoo.stab.auc)
unique(zoo.stab.auc$Taxa)

# adding mean.AUC.RR/pi and sd.AUC.RR/pi to Loop Output
#collapse between replicates here. 
zoo.stab.auc.prefin <- zoo.stab.auc %>% 
  filter(!Taxa %in% c("Nauplii", "Small cladoceran"))
unique(zoo.stab.auc.prefin$Taxa)


zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Asa", "Lake"] <- "Feresjörn"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Erk", "Lake"] <- "Erken"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Sva", "Lake"] <- "Storjärn"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Sko", "Lake"] <- "Erssjön"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Lake == "Bol", "Lake"] <- "Bolmen"

zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "F", "Treatment"] <- "Pulse"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "FS", "Treatment"] <- "Pulse & Press"
zoo.stab.auc.prefin[zoo.stab.auc.prefin$Treatment == "S", "Treatment"] <- "Press"


#write.csv2(zoo.stab.auc.prefin, 'AUCdata_biom.csv')

unique(zoo.stab.auc.prefin$Taxa)
ggplot(distinct(zoo.stab.auc.prefin, AUC.pi, AUC.RR, Replicate, AUC.RR,AUC.pi, dom, Experiment,Treatment, Taxa, Lake), aes(AUC.pi, AUC.RR,color = Taxa )) +
  geom_point(aes(shape = as.factor(Experiment),size = dom), alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.25) +                                      
  geom_hline(yintercept = 0, alpha = 0.25) +
  facet_grid(Lake~Treatment, scales = 'free') +
  ggtitle('Trophy Level: Zooplankton n = 20')+
  xlim(c(-18, 18)) +
  ylim(c(-20, 20)) + 
  labs(x = expression(AUC.~Delta~'pi'~'('~Compositional~Stability~')'), y = "Functional Stability (AUC.RR)", shape = 'Exp') +  
  theme_bw()+
  theme(legend.position = 'right')+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10, face = "bold", colour = "black", angle = 0, hjust = 0.4)) +
  theme(plot.title = element_text(size=14),
        plot.subtitle = element_text(size=13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))
#ggsave(plot = last_plot(), file = here('output/ZooAUC_reps.png'), width = 10, height = 8)


#### group only after disturbances or species for RADs####
zoo.stab.auc.prefin1 <- zoo.stab.auc.prefin %>% 
  dplyr::group_by(Taxa, Lake,Experiment, Treatment, dom) %>% # for single dots use mutate here!! 
  dplyr::mutate(mean.AUC.RR = mean(AUC.RR, na.rm = T),
                mean.AUC.pi = mean(AUC.pi, na.rm = T),
                sd.AUC.RR = sd(AUC.RR, na.rm = T),
                sd.AUC.pi = sd(AUC.pi, na.rm = T),
                se.AUC.RR = sd.AUC.RR/sqrt(n()),
                se.AUC.pi = sd.AUC.pi/sqrt(n())) %>%
  group_by(Taxa) %>%
  mutate(mean.dom = mean(con.pi))%>%
  filter(!Taxa %in% c("Nauplii", "Small cladoceran"))


#### plot with errorbars###
zoo.stab.auc.prefin1$Experiment[zoo.stab.auc.prefin1$Experiment == 1] <- 'spring'
zoo.stab.auc.prefin1$Experiment[zoo.stab.auc.prefin1$Experiment == 2] <- 'summer'

zoo.stab.auc.prefin1%>%
  distinct(mean.AUC.pi, mean.AUC.RR, se.AUC.pi, se.AUC.RR, Treatment, Taxa, Experiment,Lake,mean.dom)  %>%
ggplot(.,aes(mean.AUC.pi, mean.AUC.RR,color = Taxa )) +
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
  theme(axis.title.x = element_text(size = 16,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10, face = "bold", colour = "black", angle = 0, hjust = 0.4)) +
  theme(plot.title = element_text(size=14),
        plot.subtitle = element_text(size=13),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 6)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 12, face = 'bold'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))
ggsave(plot = last_plot(), file = here('output/Fig4.png'),width = 12, height = 8)

