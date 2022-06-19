predictive_check <- function(monitors, nSamp, n,
                             cmodel, samples, seed = 1,
                             outcome.name = "y"){
  set.seed(seed)

  pcSamples <- matrix(0, nSamp, n)
  
  for(i in 1:nSamp){
    for(monitor in monitors){
      cmodel[[monitor]] <- samples[i, monitor]
    }
    cmodel$simulate(outcome.name, includeData = TRUE)
    pcSamples[i, ] <- cmodel[[outcome.name]]
  }
  
  return(pcSamples)
}
