### Main ####

## Before running: 
# 1. Check your working directory, if not done, please set to the main folder 
getwd()

# 2. load package here
library(here)

#3. download data from Zenodo and store in the following folders:
dir.create(here('BEFD_createdData')) 
dir.create(here('SITES_Data')) 

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


