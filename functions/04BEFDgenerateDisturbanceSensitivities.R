#This function returns the species-specific disturbance-sensitivies depending on the chosen scenario

generateDisturbanceSensitivities <- function(nSpecies, 
                                   sensitivityDistribution, 
                                   equalSensitivity, 
                                   differentSensitivity,
                                   ratesOrder = 1:nSpecies,
                                   linearMinSensitivity,
                                   linearMaxSensitivity,
                                   meanNormDistr,
                                   sdNormDistr,
                                   meanLogNormDistr,
                                   sdLogNormDistr){
  if(sensitivityDistribution == 'equalSensitivities') {
    return(rep(equalSensitivity, times = nSpecies)[ratesOrder]) #give growth rates 
  }
  if(sensitivityDistribution == 'allEqualOneDifferent') {
    return(c(differentSensitivity, rep(equalSensitivity, times = nSpecies-1))[ratesOrder]) #give growth rates 
  }
  if(sensitivityDistribution == 'linearDistribution'){
    return(seq(linearMinSensitivity, linearMaxSensitivity, length.out = nSpecies)[ratesOrder])
  }
  if(sensitivityDistribution == 'normalDistribution'){
    return(rnorm(nSpecies, mean = meanNormDistr, sd = sdNormDistr)[ratesOrder])
  }
  if(sensitivityDistribution == 'logNormalDistribution'){
    return(rlnorm(nSpecies, meanlog = meanLogNormDistr, sdlog = sdLogNormDistr)[ratesOrder])
  }
}
