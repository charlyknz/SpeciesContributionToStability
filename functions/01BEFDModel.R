#This function returns the results of a full model run of the BEFD-model. The output-format is a 
#tibble. The function accepts parameters 

runBEFDModel <- function(tmax = 100, 
                         dt = 0.1, 
                         nSpecies = 5, 
                         initBiomass = 0.1, 
                         growthRatesDistribution = NA,
                         equalGrowthRate = NA,
                         minLinDistrGrowthRate = NA,
                         maxLinDistrGrowthRate = NA,
                         growthFunction = NA,
                         disturbanceIntensityPulse = NA,
                         disturbanceIntensityPress = NA,
                         sensitivityDistribution = NA,
                         equalSensitivity = NA, 
                         differentSensitivity = NA,
                         linearMinSensitivity = NA,
                         linearMaxSensitivity = NA,
                         meanNormDistr = NA,
                         sdNormDistr = NA,
                         meanLogNormDistr = NA,
                         sdLogNormDistr = NA,
                         competitionScenario = NA,
                         meanNormComp = NA, 
                         sdNormComp = NA,
                         distSensitivityOrder = 1:nSpecies,
                         capacity = 10,
                         t_pulse = 42,
                         t_press = 42,
                         numberOfRuns = 1)
{
  
  source(here('functions','02BEFDgrowthRates.R'))
  source(here('functions','03BEFDbiomassChange.R'))
  source(here('functions','04BEFDgenerateDisturbanceSensitivities.R'))
  source(here('functions','05BEFDgenerateCompetitionMatrix.R'))
  source(here('functions','06BEFDrk4Solver.R'))
  allRunsResults <- NULL
  allRunsParameters <- NULL
  
  #store fixed parameters
  fixedParameters <- tibble(
    parameterName = c('tmax', 
                      'dt', 
                      'nSpecies', 
                      'initBiomass', 
                      'growthRatesDistribution',
                      'equalGrowthRate',
                      'minLinDistrGrowthRate',
                      'maxLinDistrGrowthRate',
                      'growthFunction',
                      'disturbanceIntensityPulse',
                      'disturbanceIntensityPress',
                      'sensitivityDistribution',
                      'equalSensitivity', 
                      'differentSensitivity',
                      'linearMinSensitivity',
                      'linearMaxSensitivity',
                      'meanNormDistr',
                      'sdNormDistr',
                      'meanLogNormDistr',
                      'sdLogNormDistr',
                      'competitionScenario',
                      'meanNormComp', 
                      'sdNormComp',
                      'distSensitivityOrder',
                      'capacity',
                      't_pulse',
                      't_press',
                      'numberOfRuns'),
    parameterValue = c(tmax, 
                       dt, 
                       nSpecies, 
                       initBiomass, 
                       growthRatesDistribution,
                       equalGrowthRate,
                       minLinDistrGrowthRate,
                       maxLinDistrGrowthRate,
                       growthFunction,
                       disturbanceIntensityPulse,
                       disturbanceIntensityPress,
                       sensitivityDistribution,
                       equalSensitivity, 
                       differentSensitivity,
                       linearMinSensitivity,
                       linearMaxSensitivity,
                       meanNormDistr,
                       sdNormDistr,
                       meanLogNormDistr,
                       sdLogNormDistr,
                       competitionScenario,
                       meanNormComp, 
                       sdNormComp,
                       ifelse(sum(distSensitivityOrder == 1:nSpecies) == nSpecies,'random','different'),
                       capacity,
                       t_pulse,
                       t_press,
                       numberOfRuns),
    parameterCategory = c('time',
                          'time',
                          'initialization',
                          'initialization',
                          'growth',
                          'growth',
                          'growth',
                          'growth',
                          'growth',
                          'disturbanceSettings',
                          'disturbanceSettings',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'sensitivity',
                          'competition',
                          'competition',
                          'competition',
                          'sensitivity',
                          'growth',
                          'disturbanceSettings',
                          'disturbanceSettings',
                          'initialization'
                          )
  ) 


for (k in 1:numberOfRuns){
  
  #assign rmax-values according to the chosen distribution
  rmaxValues = growthRates(nSpecies = nSpecies,
                           growthRatesDistribution = growthRatesDistribution, 
                           equalRate = equalGrowthRate,
                           linearMin = minLinDistrGrowthRate, 
                           linearMax = maxLinDistrGrowthRate)
  
  #assign species-dependent disturbance sensitivities
  disturbanceSensitivity <- generateDisturbanceSensitivities(nSpecies = nSpecies, 
                                                             sensitivityDistribution = sensitivityDistribution, 
                                                             equalSensitivity = equalSensitivity, 
                                                             differentSensitivity = differentSensitivity,
                                                             ratesOrder = distSensitivityOrder,
                                                             linearMinSensitivity = linearMinSensitivity,
                                                             linearMaxSensitivity = linearMaxSensitivity,
                                                             meanNormDistr = meanNormDistr,
                                                             sdNormDistr = sdNormDistr,
                                                             meanLogNormDistr = meanLogNormDistr,
                                                             sdLogNormDistr = sdLogNormDistr)
  
  competitionMatrix <- generateCompetitionMatrix(competitionScenario = competitionScenario,
                                                 meanNormComp = meanNormComp,
                                                 sdNormComp = sdNormComp,
                                                 nSpecies = nSpecies)
  
  # create vectors which store the disturbance intensity settings + a control run where disturbance intensity is zero
  disturbanceIntensityPulseVec <- c(disturbanceIntensityPulse, 0)
  disturbanceIntensityPressVec <- c(disturbanceIntensityPress, 0)
  #set up second loop to run control with same parameter settings
  for (m in 1:2){    
    
    #create results-object
    modelResults <- as_tibble(matrix(nrow = tmax/dt + 1, ncol = nSpecies + 1))  %>% #row numbers for tmax-time steps + t0
      mutate_if(is.logical, as.numeric) #change columns to numeric to implement data produced by model
    
    #adjust names
    names(modelResults) <- c('timepoint', paste('species', 1:nSpecies, sep = ''))
    
    #hand over start conditions to object 
    modelResults[1,1] <- 0#first row, transform in wide format vector
    modelResults[1,1:nSpecies+1] <- t(initBiomass)#first row, transform in wide format vector

  #solve model
  for(i in 2:nrow(modelResults)){
    
    oldBiomass <- as.numeric(modelResults[i-1, -1])  # old Biomass as vector
    
    newBiomass <- oldBiomass + rk4step(dt = dt, state = oldBiomass,
                                       growthFunction = growthFunction, rmax = rmaxValues,
                                       capacity = capacity, compMatrix = competitionMatrix)
    
    #press disturbance (continuous removal of biomass after timepoint t_press)
    if(i*dt >= t_press){
      newBiomass <- newBiomass * (1 - disturbanceSensitivity * dt * disturbanceIntensityPressVec[m])
    }
    
    #Pulse disturbance (one-time removal of biomass at timepoint t_pulse)
    if(i*dt == t_pulse){
      newBiomass <-  oldBiomass* (1 - disturbanceSensitivity * disturbanceIntensityPulseVec[m])  # disturbance at timepoint
    }
    
    modelResults[i, ] <- t(c(i*dt, newBiomass)) #save disturbance values in my matrix
  }
  modelResults <- modelResults %>% 
    mutate(runNumber = k,
           runMeaning = c('withDisturbance','control')[m],
           survivingSpecies = sum(newBiomass > 1e-2))
  allRunsResults <- allRunsResults %>% 
    bind_rows(., modelResults)
  
}
  #store varying parameter values (drawn from distributions/generates within this function)
  competitionMatrixStore <- tibble(parameterName = paste('alpha',rep(1:nSpecies, each = nSpecies),
                                                         rep(1:nSpecies, times = nSpecies), sep = ''),
                                   parameterValue = as.character(c(competitionMatrix)),
                                   parameterCategory = 'competition')
  sensitivityStore <- tibble(parameterName = paste('sensitivitySpecies', 1:nSpecies, sep = ''),
                             parameterValue = as.character(disturbanceSensitivity),
                             parameterCategory = 'sensitivity')
  rmaxValuesStore <- tibble(parameterName = paste('rmaxSpecies', 1:nSpecies, sep = ''),
                            parameterValue = as.character(rmaxValues),
                            parameterCategory = 'growth')
  parameterValues <- fixedParameters %>% 
    bind_rows(.,competitionMatrixStore) %>% 
    bind_rows(.,sensitivityStore) %>% 
    bind_rows(.,rmaxValuesStore) %>% 
    mutate(runNumber = k)
  allRunsParameters <- allRunsParameters %>% 
    bind_rows(.,parameterValues) %>% 
    filter(!is.na(parameterValue))
  
  
}
  allRunsParameters <- allRunsParameters %>% 
    arrange(runNumber, parameterCategory)
return(list(simulationResults = allRunsResults,
            parameterValues = allRunsParameters))
}
