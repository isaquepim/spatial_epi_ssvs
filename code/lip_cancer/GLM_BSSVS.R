library(nimble)
library(SpatialEpi)

###### loading data ######

data(scotland)

d <- scotland$data

###### setting up the model ######

lipPriors <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  
  z0 ~ dbern(.5)
  z1 ~ dbern(.5)
  
})

lipCode <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  
  z0 ~ dbern(.5)
  z1 ~ dbern(.5)
  
  zb0 <- b0*z0
  zb1 <- b1*z1
  
  
  for (i in 1:N){
    
    log(theta[i]) <- zb0 + zb1*aff[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }

})

lipConsts <- list(N = length(d$cases))
lipData <- list(y = d$cases, aff = d$AFF, e = d$expected)
lipInits <- list(b0 = 0, b1 = 0, z0 = 0, z1 = 0)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)

lipPriorModel <- nimbleModel(code = lipPriors, name = "priors", constants = lipConsts,
                             data = lipData, inits = lipInits)

set.seed(123)


###### Compiling the model and MCMC ######

# nimbleMCMC compiles the model. If u want to customize your MCMC
# compile it yourself using:
#C.LipModel <- compileNimble(lipModel)

monitors <- c("b0","b1","z0","z1")

mcmc.out <- nimbleMCMC(code = lipCode, constants = lipConsts,
                       data = lipData, inits = lipInits,
                       nchains = 2, niter = 10000, summary = TRUE,
                       WAIC = TRUE, monitors = monitors, nburnin = 1000)

mcmc.prior <- nimbleMCMC(code = lipPriorModel, constants = lipConsts,
                         data = lipData, inits = lipInits,
                         nchains = 2, niter = 10000, summary = TRUE,
                         monitors = monitors, nburnin = 1000)

names(mcmc.out)
mcmc.out$summary

pp.plot(prior.samples = mcmc.prior$samples,
        posterior.samples = mcmc.out$samples)
