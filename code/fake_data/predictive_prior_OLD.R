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




lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  sigma ~ dinvgamma(1, 1)
  tau ~ dgamma(1, 1)
  g ~ dunif(-1,1)
  
  
  s[1:N] ~ dcar_proper(mu = mu[1:N], adj = adj[1:L], num = num[1:N],
                       tau = tau, gamma = g)
  
  for (i in 1:N){
    
    mu[i] <- mu0
    err[i] ~ dnorm(0,sigma)
    log(theta[i]) <- b0 + b1*aff[i] + s[i] + err[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
})


###### Model ######
lipConsts <- list(N = length(d$cases), L = length(LipAdj$adj),
                  adj = LipAdj$adj,
                  w = LipAdj$weights, num = LipAdj$num,
                  mu0 = 0)

lipData <- list(aff = d$AFF, e = d$expected)

lipInits <- list(b0 = 0, b1 = 0, tau = 1, sigma = 1)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)
           

cmodel  <- compileNimble(lipModel)
mcmc    <- buildMCMC(lipModel, monitors = c("b0","b1","sigma","tau","g"))
cmcmc   <- compileNimble(mcmc, project = lipModel)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)

nSamp <- nrow(samples)
n <- length(d$cases)
ppSamples <- matrix(0, nSamp, n)

set.seed(1)
for(i in 1:nSamp){
  cmodel[["b0"]] <- samples[i, "b0"]
  cmodel[["b1"]] <- samples[i, "b1"] 
  cmodel[["sigma"]] <- samples[i, "sigma"]
  cmodel[["tau"]] <- samples[i, "tau"]
  cmodel[["g"]] <- samples[i, "g"]
  cmodel$simulate("y", includeData = TRUE)
  ppSamples[i, ] <- cmodel[["y"]]
}

media <- ppSamples[1,]

for(i in 2:nSamp){
  media <- media + ppSamples[i,]
}
media <- media/nSamp
media
