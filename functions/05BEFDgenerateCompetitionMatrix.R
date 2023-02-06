generateCompetitionMatrix <- function(competitionScenario, nSpecies, meanNormComp, sdNormComp){
  if(competitionScenario == 'noCompetition'){
    compMatrix <- matrix(nrow = nSpecies, ncol = nSpecies)
    compMatrix[1:length(compMatrix)] <- 0
    diag(compMatrix) <- 1
    return(compMatrix)
  }
  if(competitionScenario == 'normalDistribution'){
    compMatrix <- matrix(nrow = nSpecies, ncol = nSpecies)
    compMatrix[1:length(compMatrix)] <- abs(rnorm(nSpecies^2, mean = meanNormComp, sd = sdNormComp))
    diag(compMatrix) <- 1
    return(compMatrix)
  }
}