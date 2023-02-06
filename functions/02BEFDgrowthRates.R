# Script giving growth rates for our model

growthRates <- function(nSpecies, growthRatesDistribution, equalRate,linearMin, linearMax){
  if(growthRatesDistribution == 'equalRates') {
    return(rep(equalRate, times = nSpecies)) #give growth rates 
  }
if(growthRatesDistribution == 'linearDistribution'){
  return(seq(linearMin, linearMax, length.out = nSpecies))
}
  }

# test the functions with two scenarios
#1. scenario
# growthRates(nSpecies = 10, growthScenario = 'equalRates', equalRate = 0.6)
#2.scenario
# growthRates(nSpecies = 10, growthScenario = 'linearDistribution', linearMin = 0.2, linearMax = 0.6)
