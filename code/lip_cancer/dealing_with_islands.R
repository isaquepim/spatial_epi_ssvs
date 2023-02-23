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
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
islands <- which(card(nb) == 0)

W.no_islands <- W[-islands, -islands]
cases.no_islands <- d$cases[-islands]
aff.no_islands <- d$AFF[-islands]
E.no_islands <- d$expected[-islands]

cases.islands <- d$cases[islands]
aff.islands <- d$AFF[islands]
E.islands <- d$expected[islands]

#matrix to adjacency list
LipAdj <- as.carAdjacency(W.no_islands)

N.no_islands <- length(cases.no_islands)
N.islands <- length(islands)

lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  tau.err ~ dgamma(3.2,1.8)
  tau ~ dgamma(1, 1)
  
  
  s[1:Ni] ~ dcar_normal(adj[1:L], w[1:L], num[1:Ni], tau, zero_mean = 1)
  
  for (i in 1:Ni){
    
    err.i[i] ~ dnorm(0,tau.err)
    log(theta.i[i]) <- log(e.i[i]) + b0 + b1*aff.i[i] + s[i] + err.i[i]
    y.i[i] ~ dpois(theta.i[i])
    
  }
  
  for (i in 1:I){
    err[i] ~ dnorm(0,tau.err)
    log(theta[i]) <- log(e[i]) + b0 + b1*aff[i] + err[i]
    y[i] ~ dpois(theta[i])
    
  }

})


###### Model ######
lipConsts <- list(Ni = N.no_islands, L = length(LipAdj$adj),
                  adj = LipAdj$adj, w = LipAdj$weights, num = LipAdj$num,
                  I = N.islands)

lipData <- list(y = cases.islands, aff = aff.islands, e = E.islands,
                y.i = cases.no_islands, aff.i = aff.no_islands, e.i = E.no_islands)

lipInits <- list(b0 = 0, b1 = 0, tau = 1, tau.err = 1, s = rep(0,N.no_islands),
                 err = rnorm(E.islands, 0,1), err.i = rnorm(N.no_islands,0,1))


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)


###### Compiling the model and MCMC ######
monitors <- c("b0","b1","tau","tau.err")

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


