#### R script to calculate the area under the curve of species-specific changes in biomass and proportion  ####
# by Charlotte Kunze 

rm(list = ls())

#load packages and functions
library(tidyverse)
library(MESS) #https://rdrr.io/cran/MESS/man/auc.html
library(ggpubr)
library(psych)
library(cowplot)
library(here)

#### import data ####
data3<- read.csv('BEFD_createdData/LRRData1.csv')
str(data3)


#### Biomass Figure ####
str(data3)
levels(as.factor(data3$species))
Lim1LRR<-ggplot(subset(data3, Limit == 'Limit1' & runNumber == 3), aes(x=timepoint, y=growth,col=species)) +
  geom_line()+
  scale_color_manual(values = c('#D56060', '#D56060', '#D56060', '#D56060', '#D56060'))+
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
  labs( x= ' ', y = "Biomass")+
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

cowplot::plot_grid(Lim1LRR, Lim2LRR, Lim3LRR,ncol = 1, nrow = 3, vjust = 1.7,hjust=0.05,label_size = 10,labels = c('a)', 'b)', 'c)'), rel_widths = c(1,1))
ggsave(plot = last_plot(), width = 8, height = 7, file = here('OutputSubmission/FigS1_Biomass.tiff'))


### create USI ####
auc <- data3 %>% 
  mutate(USI = paste(Limit, species,runNumber,Model, sep = "_")) %>% #create unique ID
  drop_na(RR)
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
names(auc)
#### calculate Area Under the Curve  ####
stab.auc <- data.frame()#opens empty data frame


# the following loop cycles through all unique case
#### AUC Loop ####
for(i in 1:length(USI)){
  temp<-auc[auc$USI==USI[i], ]#creates a temporary data frame for each case
  {
    AUC.pi.spline<- auc(temp$timepoint, temp$delta.pi, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T), type = c("spline"),absolutearea = FALSE, na.rm = T)
    AUC.RR.spline<-auc(temp$timepoint, temp$RR,  from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),type = c("spline"),absolutearea = FALSE, na.rm = FALSE)
    AUC.totRR.spline<-auc(temp$timepoint,temp$deltabm.tot, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),
                   type = c("spline"),absolutearea = FALSE,  na.rm = FALSE)
    AUC.pi<- auc(temp$timepoint, temp$delta.pi, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T), type = c("linear"),absolutearea = FALSE, na.rm = T)
    AUC.RR<-auc(temp$timepoint, temp$RR,  from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),type = c("linear"),absolutearea = FALSE, na.rm = FALSE)
    AUC.totRR<-auc(temp$timepoint,temp$deltabm.tot, from = min(temp$timepoint, na.rm = T), to = max(temp$timepoint,na.rm = T),
                   type = c("linear"),absolutearea = FALSE,  na.rm = FALSE)
    mean.con.pi <- mean(temp$con.pi)
    stab.auc<-rbind(stab.auc,
                    data.frame(temp[1,c(3,4,5,7,9,10,20)],
                               AUC.RR.spline,
                               AUC.pi.spline,
                               AUC.totRR.spline,
                               AUC.RR,
                               AUC.pi,
                               AUC.totRR,mean.con.pi))
  }
}

summary(stab.auc)

# adds sector information to dataset  
stab.auc$Sector<-1
stab.auc$Sector[stab.auc$mean.delta.pi>0&stab.auc$mean.RR>0]<-2
stab.auc$Sector[stab.auc$mean.delta.pi>0&stab.auc$mean.RR<0]<-3
stab.auc$Sector[stab.auc$mean.delta.pi<0&stab.auc$mean.RR<0]<-4


names(stab.auc)
levels(as.factor(auc$Limit))

#### save data ####
write.csv(stab.auc, 'BEFD_createdData/StabAlphaAUC.csv' )

