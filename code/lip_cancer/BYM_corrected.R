library(nimble)
library(spdep)
library(SpatialEpi)
library(sp)
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

#matrix to adjacency list
LipAdj <- as.carAdjacency(LipAdj)


lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  tau.err ~ dgamma(3.2,1.8)
  tau ~ dgamma(1, 1)
  
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau, zero_mean = 1)
  
  for (i in 1:N){
    
    err[i] ~ dnorm(0,tau.err)
    log(theta[i]) <- b0 + b1*aff[i] + s[i] + err[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
  
  sd.h <- sd(err[1:N])
  sd.c <- sd(s[1:N])
  alpha <- sd.c/(sd.h+sd.c)
})


###### Model ######
lipConsts <- list(N = length(d$cases), L = length(LipAdj$adj),
                  adj = LipAdj$adj,
                  w = LipAdj$weights, num = LipAdj$num)

lipData <- list(y = d$cases,        aff = d$AFF,
                e = d$expected)

lipInits <- list(b0 = 0, b1 = 6, tau = 1, tau.err = 1, s = rep(0,length(d$cases)),
                 err = rnorm(length(d$cases), 0,1))


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)


###### Compiling the model and MCMC ######
monitors <- c("b0","b1","tau","tau.err", "alpha")

cmodel  <- compileNimble(lipModel)
mcmc    <- buildMCMC(lipModel, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = lipModel)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)

samplesSummary(samples)

y_posterior <- predictive_check(monitors = monitors,
                                nSamp = nrow(samples),
                                n = length(d$cases),
                                cmodel = cmodel,
                                samples = samples)


d$err <- colMeans(y_posterior) - d$cases
lw <-  nb2listw(nb, style = "B", zero.policy=TRUE)
moran.test(d$err,lw, alternative="two.sided", zero.policy = TRUE)



codigo <- nimbleCode({
  b ~ dnorm(0,0.25)
})

initis <- list(b=0)

modelo <- nimbleModel(code = codigo, name = "lip", inits = initis)



###### Compiling the model and MCMC ######
monitors <- c("b")

cmodel  <- compileNimble(modelo)
mcmc    <- buildMCMC(modelo, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = modelo)
samples <- runMCMC(cmcmc, niter = 5000)
