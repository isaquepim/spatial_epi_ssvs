library(nimble)
library(sp)
library(SpatialEpi)
library(spdep)
library(ggplot2)

###### loading data ######

data(scotland)
d <- scotland$data

map <- scotland$spatial.polygon
proj4string(map) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"

map <- spTransform(map,
                   CRS("+proj=longlat +datum=WGS84 +no_defs"))

#shapefile to polygon dataframe
rownames(d) <- d$county
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)

#neighborhood list
nb <- poly2nb(map)

#nblist to matrix
LipAdj <- nb2mat(nb, style = "B", zero.policy = TRUE)
non.zero <- LipAdj[,-(which(colSums(LipAdj)==0))]
non.zero <- non.zero[-(which(rowSums(non.zero)==0)),]
non.zero.adj <- as.carAdjacency(non.zero)
#matrix to adjacency list
LipAdj <- as.carAdjacency(LipAdj)


lipPriors <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  sigma ~ dinvgamma(10, 1)
  tau ~ dgamma(1, 1)
  z0 ~ dbern(.5)
  z1 ~ dbern(.5)
  
  
})



lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  sigma ~ dinvgamma(10, 1)
  tau ~ dgamma(1, 1)
  z0 ~ dbern(.5)
  z1 ~ dbern(.5)
  
  zb0 <- b0*z0
  zb1 <- b1*z1
  
  # stand-alone priors for comparison
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau, zero_mean = 1)
  
  for (i in 1:N){
    
    err[i] ~ dnorm(0,sigma)
    log(theta[i]) <- zb0 + zb1*aff[i] + s[i] + err[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
})


###### Model ######
lipConsts <- list(N = length(d$cases), L = length(LipAdj$adj),
                  adj = LipAdj$adj,
                  w = LipAdj$weights, num = LipAdj$num)

lipData <- list(y = d$cases,        aff = d$AFF,
                e = d$expected)

lipInits <- list(b0 = 0, b1 = 0, tau = 1, sigma = 1, z0 = 1, z1 = 1)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)

lipPriorModel <- nimbleModel(code = lipPriors, name = "priors", constants = lipConsts,
                             data = lipData, inits = lipInits)

###### Compiling the model and MCMC ######

monitors <- c("b0","b1","tau",
              "sigma","z0","z1")
mcmc.out <- nimbleMCMC(code = lipCode, constants = lipConsts,
                       data = lipData, inits = lipInits,
                       nchains = 2, niter = 10000, summary = TRUE,
                       WAIC = TRUE, monitors = monitors, nburnin = 1000)
mcmc.prior <- nimbleMCMC(code = lipPriorModel, constants = lipConsts,
                         data = lipData, inits = lipInits,
                         nchains = 2, niter = 10000, summary = TRUE,
                         monitors = monitors, nburnin = 1000)

mcmc.out$summary


pp.plot(prior.samples = mcmc.prior$samples,
        posterior.samples = mcmc.out$samples)

