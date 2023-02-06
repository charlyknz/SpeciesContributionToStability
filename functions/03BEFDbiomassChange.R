# Script to set growth functions for our model
biomassChange <- function(growthFunction, biomass, rmax, capacity, compMatrix){
  if(growthFunction == 'densityDependent'){
    return(rmax*biomass*(1-biomass/capacity))
  }
  if(growthFunction == 'competitionDependent'){
    return(rmax*biomass*(1 - colSums(compMatrix * biomass  / capacity)))
  }
}