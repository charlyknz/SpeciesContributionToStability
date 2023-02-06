# Model implementation species specific contribution to stability
# in this script we explore the effects of species competition on system dynamics
# by Dominik Bahlburg 
setwd("~/Desktop/BEFD22/BEFDisturbance/Untitled/")

#load packages and functions
library(tidyverse)
library(here)
source(here('functions/01BEFDModel.R'))

#------------------------------------------------------------------------#
# A few information before you run the model:
# In general, all model parameters are adjustable. However, for simplicity we assume that some of the
# parameters are fixed (e.g. number of timesteps). 
# These parameters have default values and are not changed if not not explicitly 
# stated otherwise.
# Changes can be made every time the model is executed
#
# These fixed parameters and their default values are listed below:
#
# capacity: 10, biomass concentration at capacity
# nSpecies: number of species to simulate

# Time variables
# tmax: 500, stimulation steps - timepoints which will be observed
# dt: 0.1, step size, Number of simulated timepoints: tmax/dt
# t_pulse: 42, Timepoint when pulse disturbance happens
# t_press: 42, Timepoint when press disturbance happens
# numberOfRuns: 1, how often should the model run? Important for runs which include stochasticity

# Biomass at timepoint t0
# initBiomass: 0.1, biomass concentration for each species at t0
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
# arguments/parameters that that need to be passed to the function at each model run:
# When scenarios are chosen such as "equalRates", all additional scenario-dependent
# parameters are listed below as "+parameter"

# growthRatesDistribution. Options: "equalRates", "linearDistribution"
# in case "equalRates": +equalGrowthRate, 
# in case "linearDistribution": +minLinDistrGrowthRate +maxLinDistrGrowthRate
#
# growthFunction. Options: "densityDependent", "competitionDependent"
# in case "competitionDependent: +competitionScenario
#         in case "normalDistribution": +meanNormComp, +sdNormComp
#
# disturbanceIntensity: numerical value between 0 and 1 for either pulse, press or pulse and press in combination
# disturbanceIntensityPulse = 0.5
# disturbanceIntensityPress = 0.0025
#
# sensitivityDistribution. Options: "equalSensitivities", "allEqualOneDifferent" , "linearDistribution", "normalDistribution", "logNormalDistribution"
# in case "equalSensitivities": +equalSensitivity
# in case "allEqualOneDifferent": +equalSensitivity, +differentSensitivity
# in case "linearDistribution": +linearMinSensitivity, +linearMaxSensitivity
# in case "normalDistribution": +meanNormDistr, +sdNormDistr
# in case "logNormalDistribution": +meanLogNormDistr,+sdLogNormDistr
# 
#
#Limit scenarios:
# Limit 1: All species respond equally to the disturbance 
# Limit 2: All species are equally sensitive to the disturbance but one is resistant/ benefits
# Limit 3: All species are equally resistant to the disturbance but one is sensitive/ struggles
#------------------------------------------------------------------------#
# Run/test the model with test

#### Limit 1 ####
allSensitive_pulse <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'equalSensitivities',
                                       equalSensitivity = 1,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0.5,
                                       disturbanceIntensityPress = 0.0,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
allSensitive_pulseResults <- allSensitive_pulse[[1]]
parametersallSensitive_pulse <- allSensitive_pulse[[2]]

allSensitive_pulseResults <- allSensitive_pulseResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(allSensitive_pulse, file=here('modelOutput','Limit1ScenarioallSensitive_pulse.RData'))

#Load List object
load(here('modelOutput','Limit1ScenarioallSensitive_pulse.RData'))

#export results
#write.csv(allSensitive_pulseResults, here('modelOutput','Limit1SensitiveResults.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- allSensitive_pulse[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit1All <- left_join(allSensitive_pulseResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit1'),
         Model = 'pulse')
#Save simulation results as csv
write.csv(Limit1All, here('modelOutput/Loop','Limit1AlphaModelResults_pulse.csv'))


### pulse and press ###
allSensitive_pulsepresspress <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                                 sensitivityDistribution = 'equalSensitivities',
                                                 equalSensitivity = 1,
                                                 growthFunction = 'competitionDependent',
                                                 equalGrowthRate = 0.6,
                                                 disturbanceIntensityPulse = 0.5,
                                                 disturbanceIntensityPress = 0.0025,
                                                 t_pulse = 150,
                                                 t_press = 150,
                                                 competitionScenario = 'normalDistribution',
                                                 meanNormComp = 0.7, 
                                                 sdNormComp = 0.05,
                                                 tmax = 450,
                                                 nSpecies = 5,
                                                 numberOfRuns = 50,
                                                 dt = 0.5)

#extract simulation results and parameter values
allSensitive_pulsepresspressResults <- allSensitive_pulsepresspress[[1]]
parametersallSensitive_pulsepresspress <- allSensitive_pulsepresspress[[2]]

allSensitive_pulsepresspressResults <- allSensitive_pulsepresspressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)


#Save list object
save(allSensitive_pulsepresspress, file=here('modelOutput','Limit1ScenarioallSensitive_pulsepresspress.RData'))

#Load List object
load(here('modelOutput','Limit1ScenarioallSensitive_pulsepresspress.RData'))

#export results
#write.csv(allSensitive_pulsepresspressResults, here('modelOutput','Limit1SensitiveResults_pulsepress.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- allSensitive_pulsepresspress[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit1All <- left_join(allSensitive_pulsepresspressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit1'),
         Model = 'pulsepress')
#Save simulation results as csv
write.csv(Limit1All, here('modelOutput/Loop','Limit1AlphaModelResults_pulsepress.csv'))


### press ###
allSensitive_press <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'equalSensitivities',
                                       equalSensitivity = 1,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0,
                                       disturbanceIntensityPress = 0.0025,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
allSensitive_pressResults <- allSensitive_press[[1]]
parametersallSensitive_press <- allSensitive_press[[2]]

allSensitive_pressResults <- allSensitive_pressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(allSensitive_press, file=here('modelOutput','Limit1ScenarioallSensitive_press.RData'))

#Load List object
load(here('modelOutput','Limit1ScenarioallSensitive_press.RData'))

#export results
#write.csv(allSensitive_pressResults, here('modelOutput','Limit1SensitiveResults_press.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- allSensitive_press[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit1All <- left_join(allSensitive_pressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit1'),
         Model = 'press')
#Save simulation results as csv
write.csv(Limit1All, here('modelOutput/Loop','Limit1AlphaModelResults_press.csv'))





#### Limit 2 ####
oneLessSensitive_pulse <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                 sensitivityDistribution = 'allEqualOneDifferent',
                                 equalSensitivity = 1,
                                 differentSensitivity = 0.5,
                                 growthFunction = 'competitionDependent',
                                 equalGrowthRate = 0.6,
                                 disturbanceIntensityPulse = 0.5,
                                 disturbanceIntensityPress = 0.0,
                                 t_pulse = 150,
                                 t_press = 150,
                                 competitionScenario = 'normalDistribution',
                                 meanNormComp = 0.7, 
                                 sdNormComp = 0.05,
                                 tmax = 450,
                                 nSpecies = 5,
                                 numberOfRuns = 50,
                                 dt = 0.5)

#extract simulation results and parameter values
oneLessSensitive_pulseResults <- oneLessSensitive_pulse[[1]]
parametersoneLessSensitive_pulse <- oneLessSensitive_pulse[[2]]

oneLessSensitive_pulseResults <- oneLessSensitive_pulseResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(oneLessSensitive_pulse, file=here('modelOutput','Limit2ScenarioOneLessSensitive_pulse.RData'))

#Load List object
load(here('modelOutput','Limit2ScenarioOneLessSensitive_pulse.RData'))

#export results
#write.csv(oneLessSensitive_pulseResults, here('modelOutput','Limit2SensitiveResults.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneLessSensitive_pulse[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit2All <- left_join(oneLessSensitive_pulseResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit2'),
  Model = 'pulse')
#Save simulation results as csv
write.csv(Limit2All, here('modelOutput/Loop','Limit2AlphaModelResults_pulse.csv'))


### pulse and press ###
oneLessSensitive_pulsepresspress <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'allEqualOneDifferent',
                                       equalSensitivity = 1,
                                       differentSensitivity = 0.5,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0.5,
                                       disturbanceIntensityPress = 0.0025,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
oneLessSensitive_pulsepresspressResults <- oneLessSensitive_pulsepresspress[[1]]
parametersoneLessSensitive_pulsepresspress <- oneLessSensitive_pulsepresspress[[2]]

oneLessSensitive_pulsepresspressResults <- oneLessSensitive_pulsepresspressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)


#Save list object
save(oneLessSensitive_pulsepresspress, file=here('modelOutput','Limit2ScenarioOneLessSensitive_pulsepresspress.RData'))

#Load List object
load(here('modelOutput','Limit2ScenarioOneLessSensitive_pulsepresspress.RData'))

#export results
#write.csv(oneLessSensitive_pulsepresspressResults, here('modelOutput','Limit2SensitiveResults_pulsepress.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneLessSensitive_pulsepresspress[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit2All <- left_join(oneLessSensitive_pulsepresspressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit2'),
         Model = 'pulsepress')
#Save simulation results as csv
write.csv(Limit2All, here('modelOutput/Loop','Limit2AlphaModelResults_pulsepress.csv'))


### press ###
oneLessSensitive_press <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'allEqualOneDifferent',
                                       equalSensitivity = 1,
                                       differentSensitivity = 0.5,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0,
                                       disturbanceIntensityPress = 0.0025,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
oneLessSensitive_pressResults <- oneLessSensitive_press[[1]]
parametersoneLessSensitive_press <- oneLessSensitive_press[[2]]

oneLessSensitive_pressResults <- oneLessSensitive_pressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(oneLessSensitive_press, file=here('modelOutput','Limit2ScenarioOneLessSensitive_press.RData'))

#Load List object
load(here('modelOutput','Limit2ScenarioOneLessSensitive_press.RData'))

#export results
#write.csv(oneLessSensitive_pressResults, here('modelOutput','Limit2SensitiveResults_press.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneLessSensitive_press[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit2All <- left_join(oneLessSensitive_pressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit2'),
         Model = 'press')
#Save simulation results as csv
write.csv(Limit2All, here('modelOutput/Loop','Limit2AlphaModelResults_press.csv'))





#### Limit 3 ####

# Run/test the model with test
oneMoreSensitive_pulse <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'allEqualOneDifferent',
                                       equalSensitivity = 0.5,
                                       differentSensitivity = 1,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0.5,
                                       disturbanceIntensityPress = 0.0,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
oneMoreSensitive_pulseResults <- oneMoreSensitive_pulse[[1]]
parametersoneMoreSensitive_pulse <- oneMoreSensitive_pulse[[2]]

oneMoreSensitive_pulseResults <- oneMoreSensitive_pulseResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(oneMoreSensitive_pulse, file=here('modelOutput','Limit3ScenariooneMoreSensitive_pulse.RData'))

#Load List object
load(here('modelOutput','Limit3ScenariooneMoreSensitive_pulse.RData'))

#export results
#write.csv(oneMoreSensitive_pulseResults, here('modelOutput','Limit3SensitiveResults.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneMoreSensitive_pulse[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit3All <- left_join(oneMoreSensitive_pulseResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit3'),
         Model = 'pulse')
#Save simulation results as csv
write.csv(Limit3All, here('modelOutput/Loop','Limit3AlphaModelResults_pulse.csv'))


### pulse and press ###
oneMoreSensitive_pulsepresspress <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                                 sensitivityDistribution = 'allEqualOneDifferent',
                                                 equalSensitivity = 0.5,
                                                 differentSensitivity = 1,
                                                 growthFunction = 'competitionDependent',
                                                 equalGrowthRate = 0.6,
                                                 disturbanceIntensityPulse = 0.5,
                                                 disturbanceIntensityPress = 0.0025,
                                                 t_pulse = 150,
                                                 t_press = 150,
                                                 competitionScenario = 'normalDistribution',
                                                 meanNormComp = 0.7, 
                                                 sdNormComp = 0.05,
                                                 tmax = 450,
                                                 nSpecies = 5,
                                                 numberOfRuns = 50,
                                                 dt = 0.5)

#extract simulation results and parameter values
oneMoreSensitive_pulsepresspressResults <- oneMoreSensitive_pulsepresspress[[1]]
parametersoneMoreSensitive_pulsepresspress <- oneMoreSensitive_pulsepresspress[[2]]

oneMoreSensitive_pulsepresspressResults <- oneMoreSensitive_pulsepresspressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)


#Save list object
save(oneMoreSensitive_pulsepresspress, file=here('modelOutput','Limit3ScenariooneMoreSensitive_pulsepresspress.RData'))

#Load List object
load(here('modelOutput','Limit3ScenariooneMoreSensitive_pulsepresspress.RData'))

#export results
#write.csv(oneMoreSensitive_pulsepresspressResults, here('modelOutput','Limit3SensitiveResults_pulsepress.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneMoreSensitive_pulsepresspress[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit3All <- left_join(oneMoreSensitive_pulsepresspressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit3'),
         Model = 'pulsepress')
#Save simulation results as csv
write.csv(Limit3All, here('modelOutput/Loop','Limit3AlphaModelResults_pulsepress.csv'))


### press ###
oneMoreSensitive_press <- runBEFDModel(growthRatesDistribution = 'equalRates',
                                       sensitivityDistribution = 'allEqualOneDifferent',
                                       equalSensitivity = 0.5,
                                       differentSensitivity = 1,
                                       growthFunction = 'competitionDependent',
                                       equalGrowthRate = 0.6,
                                       disturbanceIntensityPulse = 0,
                                       disturbanceIntensityPress = 0.0025,
                                       t_pulse = 150,
                                       t_press = 150,
                                       competitionScenario = 'normalDistribution',
                                       meanNormComp = 0.7, 
                                       sdNormComp = 0.05,
                                       tmax = 450,
                                       nSpecies = 5,
                                       numberOfRuns = 50,
                                       dt = 0.5)

#extract simulation results and parameter values
oneMoreSensitive_pressResults <- oneMoreSensitive_press[[1]]
parametersoneMoreSensitive_press <- oneMoreSensitive_press[[2]]

oneMoreSensitive_pressResults <- oneMoreSensitive_pressResults%>%
  gather( key = 'species', value = 'growth', -timepoint, -runNumber,-runMeaning, -survivingSpecies)

#Save list object
save(oneMoreSensitive_press, file=here('modelOutput','Limit3ScenariooneMoreSensitive_press.RData'))

#Load List object
load(here('modelOutput','Limit3ScenariooneMoreSensitive_press.RData'))

#export results
#write.csv(oneMoreSensitive_pressResults, here('modelOutput','Limit3SensitiveResults_press.csv'))

#Extract alpha- and sensitivity-values for post-simulation analysis
#first: extract list element "parameterValues" and filter competition- and sensitivity-related parameters
#create new column which contains the species name and, in case of alpha, the competitor
parameterValues <- oneMoreSensitive_press[['parameterValues']] %>%
  filter(parameterCategory %in% c('competition','sensitivity')) %>%
  filter(str_detect(parameterName, 'alpha|sensitivitySpecies')) %>%
  mutate(species = paste('species',ifelse(str_detect(parameterName, 'alpha'),
                                          substr(parameterName, 6, 6),
                                          substr(parameterName, 19, 19)),sep = ''),
         competitor =  ifelse(str_detect(parameterName, 'alpha'),
                              paste('species',substr(parameterName, 7, 7),sep = ''),
                              NA),
         parameterValue = as.numeric(parameterValue))

#summarise alpha values -> mean and sd for alpha values across species excluding(!)
#intraspecific competition (removed with filter)
alphaSummary <- parameterValues %>%
  # filter(species != competitor) %>%
  group_by(parameterCategory, runNumber, species) %>%
  summarise(meanValue = mean(parameterValue)) %>%
  spread(parameterCategory, meanValue)


#combine alpha and sensitivity with model output 
Limit3All <- left_join(oneMoreSensitive_pressResults, alphaSummary, by = c('species', 'runNumber'))%>%
  mutate(Limit = paste('Limit3'),
         Model = 'press')
#Save simulation results as csv
write.csv(Limit3All, here('modelOutput/Loop','Limit3AlphaModelResults_press.csv'))





#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#### loop importing single model runs ####
#Rausziehen der Dateipfade:
setwd("~/Desktop/BEFD22/BEFDisturbance/Untitled/modelOutput/Loop")
directory <- paste(getwd(),paste= '',sep = '') # hier sind meine Daten
files <- dir(directory, recursive=TRUE, full.names=TRUE, pattern="\\.csv$")
liste = list.files(pattern = "*.csv")

#Erstellen eines leeren "Blanko-Objektes", welches spaeter alle Daten speichern soll
data <- NULL
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#
directory <- paste(getwd(),'/modelOutput/Loop',sep = ';') # hier sind meine Daten
liste = list.files(pattern = "*.csv")


tibble <- list.files(pattern = "*.csv") %>%  #lists all files with csv ending in our wd
  map_df(~read_csv(.)) #%>% #import data
  #filter(survivingSpecies == 5)
#separate control data to rename columns
all.con <-tibble%>%
  filter(runMeaning != 'withDisturbance') %>%
  rename(con.bio = growth) %>%
  select(-runMeaning, -...1)
names(all.con)
limit_sameR<- tibble %>%
  filter(runMeaning == 'withDisturbance') %>%
  select(-runMeaning, -...1)


names(limit_sameR)
LimitDataConDist <- right_join(limit_sameR, all.con, by = c('timepoint', 'runNumber', 'survivingSpecies', 'species', 'Limit','competition', 'sensitivity', 'Model')) 

# for comparison, calculate the overall community response
data.tot<-LimitDataConDist %>%
  group_by(timepoint, Limit, runNumber, Model) %>% 
  summarise(treat.tot = sum(growth), #sum of species growth 
            con.tot = sum(con.bio))

summary(data.tot)

# LRR of community 
data.tot$LRR.tot<-log(data.tot$treat.tot/data.tot$con.tot)
data.tot$deltabm.tot<-(data.tot$treat.tot-data.tot$con.tot)/(data.tot$treat.tot+data.tot$con.tot)
names(data.tot)
hist(data.tot$LRR.tot)
hist(data.tot$deltabm.tot)
str(data.tot )

#create data set for species specific effect sizes
names(LimitDataConDist)
data3<-NULL
data3<-LimitDataConDist %>%
  left_join(.,data.tot, by = c("timepoint","Limit", 'runNumber', 'Model'))
summary(data3)

#take out those rows where biomass is 0 in both treatment + control
data3$RR<-data3$growth+data3$con.bio
data3<-filter(data3, RR!=0)

# create species specific LRR for biomass 
data3$LRR<-log(data3$growth/data3$con.bio)
#create species specific change in biomass
data3$RR<-(data3$growth-data3$con.bio)/(data3$growth+data3$con.bio)
#the absence of species in control or treatment creates INF, get rid of these
data3$LRR[data3$LRR=="Inf"]<-NA
data3$LRR[data3$LRR=="-Inf"]<-NA

# create species specific contribution to biomass in each treatment
data3$treat.pi<-data3$growth/data3$treat.tot
data3$con.pi<-data3$con.bio/data3$con.tot

# create species specific LRR and difference for pi
data3$delta.pi<-data3$treat.pi-data3$con.pi
which(is.na(data3$delta.pi))

#weight LRR  by mean pi
data3$mean.pi<-0.5*(data3$treat.pi+data3$con.pi)
data3$LRR.w<-data3$LRR*data3$mean.pi
data3$RR.w<-data3$RR*data3$mean.pi

#create the deviance between species and community effect sizes
data3$LRR.diff<-data3$LRR-data3$LRR.tot
data3$LRR.diff.w<-data3$LRR.diff*data3$mean.pi
summary(data3)
str(data3)

#### create First graphs ####
#plot
delta.pi.plot<-ggplot(data3, aes(x=timepoint, y=delta.pi,
                                 col=species)) +
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  theme(legend.position="bottom")+
  # scale_y_continuous(limit = c(0.15, -0.1), breaks = seq(-0.1,0.1,0.05))  +
  geom_line()+
  xlab("time") +  ylab ("Delta.pi")+
  facet_grid(~Limit~Model)+
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) 
delta.pi.plot
#ggsave(plot =delta.pi.plot, file = 'sameRnoCompSpecies_contribution_mean.jpeg', width=8, height = 5)


LRR.plot<-ggplot(data3, aes(x=timepoint, y=LRR,
                            col=species, alpha = timepoint)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept=0, col="grey")+
  geom_vline(xintercept=0, col="grey")+
  theme(legend.position="bottom")+
  # scale_y_continuous(limit = c(0.15, -0.1), breaks = seq(-0.1,0.1,0.05))  +
  geom_line()+
  xlab("time") +  ylab ("LRR biomass")+ ggtitle('LRR of multiple runs')+
  facet_grid(~Model)
LRR.plot
ggsave(plot =LRR.plot, file = 'sameRnoCompSpeciesLRR_mean.jpeg', width=8, height = 5)

#### write csv ####
write.csv2(data3,'LRRData2.csv')

