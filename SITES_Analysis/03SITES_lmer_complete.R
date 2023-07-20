### for Charly

rm(list=ls())
graphics.off()

library(tidyverse)
library(lme4)
library(sjPlot)

#### read data ####
data<-read.csv2("~/Desktop/phD/SpeciesContributionToStability/complete/AUCdata_3.csv")
summary(data)
unique(data$Taxa)

#### Model formulation for each disturbance type and dimension ####
# we introduced Lake and Experiment as random effects

mod1.press.aucRR<-lmer(AUC.RR~Taxa-1+ (1|Lake) + (1|Experiment),
                       data=data[data$Treatment=="Press",])
summary(mod1.press.aucRR)

mod2.press.aucPI<-lmer(AUC.pi~Taxa-1+ (1|Lake) + (1|Experiment) ,
                       data=data[data$Treatment=="Press",])
summary(mod2.press.aucPI)


mod3.pulse.aucRR<-lmer(AUC.RR~Taxa-1 + (1|Lake) + (1|Experiment),
                       data=data[data$Treatment=="Pulse",])
summary(mod3.pulse.aucRR)

mod4.pulse.aucPI<-lmer(AUC.pi~Taxa-1 + (1|Lake) + (1|Experiment),
                       data=data[data$Treatment=="Pulse",])
summary(mod4.pulse.aucPI)


mod5.presspulse.aucRR<-lmer(AUC.RR~Taxa-1 + (1|Lake) + (1|Experiment),
                       data=data[data$Treatment=="Pulse & Press",])
summary(mod5.presspulse.aucRR)

mod6.presspulse.aucPI<-lmer(AUC.pi~Taxa-1 + (1|Lake) + (1|Experiment),
                       data=data[data$Treatment=="Pulse & Press",])
summary(mod6.presspulse.aucPI)


#### summary table for the model outputs ####
tab_model(mod1.press.aucRR,mod2.press.aucPI,mod3.pulse.aucRR,
          mod4.pulse.aucPI,mod5.presspulse.aucRR,
          mod6.presspulse.aucPI,digits = 4, show.ci = FALSE, file = here("complete/YOURTABLENAME.doc"))
