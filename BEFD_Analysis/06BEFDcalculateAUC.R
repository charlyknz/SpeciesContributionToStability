#### R script to calculate AUC  ####
# by Charlotte Kunze 

#load packages and functions
library(reshape2)
library(tidyverse)
library(MESS) #https://rdrr.io/cran/MESS/man/auc.html
library(nlme)
library(ggpubr)
library(psych)
library(cowplot)


#### import data ####
data3<- read.csv2('BEFD_createdData/LRRData2.csv')
str(data3)


#### Fig. 1####
str(data3)
levels(as.factor(data3$species))
Lim1LRR<-ggplot(subset(data3, Limit == 'Limit1' & runNumber == 3), aes(x=timepoint, y=growth,col=species)) +
  geom_line()+
  scale_color_manual(values = c('#68789E', '#68789E', '#68789E', '#68789E', '#68789E'))+
  labs(x= ' ',  y = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 12)) +
  theme(legend.position="none")
Lim1LRR

Lim2LRR<-ggplot(subset(data3, Limit == 'Limit2' & runNumber == 40), aes(x=timepoint, y=growth,
                                                                  col=species)) +
  scale_color_manual(values = c('#68789E', '#C6C6C680', '#C6C6C680', '#C6C6C680', '#C6C6C680'))+
  geom_line()+
  labs( x= ' ', y = "biomass")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 12)) +
  theme(legend.position="none")
Lim2LRR

Lim3LRR<-ggplot(subset(data3, Limit == 'Limit3'&runNumber == 27), aes(x=timepoint, y=growth,
                                                                  col=species)) +
  geom_line()+
  scale_color_manual(values = c('#D56060', '#C6C6C680', '#C6C6C680', '#C6C6C680', '#C6C6C680'))+
  labs(x = 'time', y = " ")+
  facet_wrap(~Model, ncol = 3)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 12)) +  theme(legend.position="none")
Lim3LRR

plot_grid(Lim1LRR, Lim2LRR, Lim3LRR,ncol = 1, nrow = 3, vjust = 1.7,hjust=0.05,label_size = 10,labels = c('a)', 'b)', 'c)'), rel_widths = c(1,1))
ggsave(plot = last_plot(), width = 8, height = 7, file = here('output/Fig.1_Biomass.png'))


### create USI ####
auc <- data3 %>% 
  mutate(USI = paste(Limit, species,runNumber,Model, sep = "_")) %>% #create unique ID
  drop_na(LRR)
which(is.infinite(auc$deltabm.tot))
str(auc)

#order data after time
auc <-auc[with(auc, order(USI, timepoint)),]

#First sampling needs to be 0 in order to use the AUC function
#auc$timepoint <- ifelse(auc$timepoint == 0.0, 0, auc$timepoint)
str(auc)
#create vector with the unique ID
USI<-unique(auc$USI)# USI is the identifier for the experimental unit 

head(USI)

#### calculate Area Under the Curve  ####
stab.auc <- data.frame()#opens empty data frame


# the following loop cycles through all unique case
#### AUC Loop ####
for(i in 1:length(USI)){
  temp<-auc[auc$USI==USI[i], ]#creates a temporary data frame for each case
  {
    AUC.pi<- auc(temp$timepoint, temp$delta.pi, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T), type = c("spline"),absolutearea = FALSE, na.rm = T)
    AUC.RR<-auc(temp$timepoint, temp$RR,  from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),type = c("spline"),absolutearea = FALSE, na.rm = FALSE)
    AUC.totRR<-auc(temp$timepoint,temp$deltabm.tot, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),
                   type = c("spline"),absolutearea = FALSE,  na.rm = FALSE)
    mean.delta.pi<-mean(temp$delta.pi)
    mean.RR<-mean(temp$RR)
    mean.pi<-mean(temp$mean.pi)
    mean.con.pi <- mean(temp$con.pi)
    mean.treat.pi <- mean(temp$treat.pi)
    sd.delta.pi<-sd(temp$delta.pi)
    sd.RR<-sd(temp$RR)
    stab.auc<-rbind(stab.auc,
                    data.frame(temp,
                               AUC.RR,
                               AUC.pi,
                               AUC.totRR,mean.con.pi, mean.treat.pi,
                               mean.delta.pi,mean.RR,mean.pi,
                               sd.delta.pi,sd.RR))
  }
}

summary(stab.auc)

stab.auc$Sector<-1
stab.auc$Sector[stab.auc$mean.delta.pi>0&stab.auc$mean.RR>0]<-2
stab.auc$Sector[stab.auc$mean.delta.pi>0&stab.auc$mean.RR<0]<-3
stab.auc$Sector[stab.auc$mean.delta.pi<0&stab.auc$mean.RR<0]<-4

stab.auc$RR2<-stab.auc$AUC.RR
stab.auc$RR2[stab.auc$Sector==3]<-(-stab.auc$AUC.RR[stab.auc$Sector==3])
stab.auc$RR2[stab.auc$Sector==4]<-(-stab.auc$AUC.RR[stab.auc$Sector==4])
stab.auc$pi2<-stab.auc$AUC.pi
stab.auc$pi2[stab.auc$Sector==4]<-(-stab.auc$AUC.pi[stab.auc$Sector==4])
stab.auc$pi2[stab.auc$Sector==1]<-(-stab.auc$AUC.pi[stab.auc$Sector==1])


names(stab.auc)

### calculate relative alpha ###
max <- stab.auc %>%
  mutate(competition, max = max(competition))%>%
  ungroup() %>%
  distinct( Model, runNumber,  species,max, competition, sensitivity,  mean.treat.pi, mean.con.pi, AUC.totRR,AUC.pi, AUC.RR, Limit, mean.pi,LRR, timepoint)

levels(as.factor(auc$Limit))
names(max)
#### write csv ####
write.csv(max, 'StabAlphaAUC2.csv' )

