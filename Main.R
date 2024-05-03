### Main ####

## Before running: 
# 1. Check your working directory, if not done, please set to the main folder 
getwd()


# 2. load R packages and environment used to develop the code including package versions 

install.packages('renv')
library(renv)
renv::restore()

library(here)

## Sometimes the renv function does not work (depending on the R version on your local device),
# ALTERNATIVE: install packages one by one

#install.packages('cowplot')
#install.packages('here')
#library(here)
#install.packages('ggpubr') 
#install.packages('ggpmisc') 
#install.packages('lme4')
#install.packages('MESS')
#install.packages('psych')
#install.packages('sjPlot')
#install.packages('tidyverse')


#3. Download data from Zenodo and store in the following folders:
dir.create(here('BEFD_createdData')) 
dir.create(here('SITES_Data')) 

# Data Link: https://doi.org/10.5281/zenodo.11046700

# Before running create submission folder
dir.create(here('OutputSubmission')) 

#4. Model simulations

#4.1 If you want to create model simulations yourself run the following - this takes a while!!!
#source(here("BEFD_Analysis/05BEFDcreateData.R"))

#4.2 Create simulation plots and start simulation analysis

##calculation of AUC from LRR csv file
source(here("BEFD_Analysis/06BEFDcalculateAUC.R"))

##Fig 2
source(here("BEFD_Analysis/07BEFDanalyseAUC.R"))

##Fig 3 & correlation plots for statistical analysis
source(here("BEFD_Analysis/08BEFDdominance.R"))


#5. Empirical Data

#5.1 calculation of the area und the curve (AUC) of species specific responses given by zooplankt.csv
# Creates plot Fig. 4
source(here("SITES_Analysis/01SITES_createTidyData_complete.R"))

#5.2 Creates plot Fig. 5
source(here("SITES_Analysis/02SITES_TidyanalyseTraits_complete.R"))

#5.3 Perform LMER to calculate likelihood of species contributions  
source(here("SITES_Analysis/03SITES_lmer_complete.R"))


