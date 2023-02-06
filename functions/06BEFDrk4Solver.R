rk4step <- function(dt, state, growthFunction, rmax, capacity, compMatrix){
  
  k1 <- biomassChange(biomass = state, growthFunction = growthFunction, rmax = rmax, capacity = capacity, compMatrix = compMatrix)
  k2 <- biomassChange(biomass = state+0.5*k1*dt, growthFunction = growthFunction, rmax = rmax, capacity = capacity, compMatrix = compMatrix)
  k3 <- biomassChange(biomass = state+0.5*k2*dt, growthFunction = growthFunction, rmax = rmax, capacity = capacity, compMatrix = compMatrix)
  k4 <- biomassChange(biomass = state+k3*dt, growthFunction = growthFunction, rmax = rmax, capacity = capacity, compMatrix = compMatrix)

  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}
