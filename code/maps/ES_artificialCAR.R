library(ggplot2)
library(ggmap)
library(rgdal)
library(maps)
library(mapdata)
library(broom)
library(MASS)
library(patchwork)
library(sf)
library(sp)
library(spdep)
library(rwc)
library(dplyr)

set.seed(65465)

base_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\code\\maps"


ES.shp <- readOGR(
  dsn = paste0(base_dir, '\\shapefiles'),
  layer = 'ES_G'
)

ES.fortified <- tidy(ES.shp)
ES.sf <- st_as_sf(ES.shp)


nb <- poly2nb(ES.shp)
adj <- nb2mat(nb,   style = "B", zero.policy = TRUE)
lw <-  nb2listw(nb, style = "B", zero.policy=TRUE)

betas <- c(1,-1,0.4)

x1 <- rnorm(n=78, 0,1)
x2 <- rnorm(n=78, 0,1)


ndi <- rep(0,78)

for(i in 1:78){
  ndi[i] <- sum(adj[i,])
}

rho <- 0.99
Precisao <- diag(ndi) - rho*adj
CAR <- rnorm.Q(Precisao, zero.constraint = TRUE)
ES.sf$CAR <- as.vector(CAR)

ggplot()+
  geom_sf(data = ES.sf, aes(fill = G))


moran(ES.sf$G, lw, length(nb), Szero(lw))
moran.test(ES.sf$G,lw, alternative="two.sided")

beta0 <- 1.5
beta1 <- 0.5
beta2 <- -0.5
beta3 <- 0.3

x1 <- rnorm(n=78, 0,1)
x2 <- rnorm(n=78, 0,1)

ES.sf$x1 <- x1
ES.sf$x2 <- x2

log.mu <- beta0 + ES.sf$x1*beta1 + ES.sf$x2*beta2 + ES.sf$G*beta3 
mu <- exp(log.mu)
mu
y <- rpois(n=78,mu)

ES.sf$Y <- y

ggplot()+
  geom_sf(data = ES.sf, aes(fill = Y))




library(nimble)

ESCode <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  b2 ~ dnorm(0,1)
  
  for (i in 1:N){
    
    log(theta[i]) <- b0 + b1*x1[i] + b2*x2[i]
    y[i] ~ dpois(theta[i])
    
  }
})


ESConsts <- list(N = 78)
ESData <- list(y = ES.sf$Y, x1 = ES.sf$x1, x2 = ES.sf$x2)
ESInits <- list(b0 = 0, b1 = 0, b2 = 0)


ESModel <- nimbleModel(code = ESCode, name = "ES", constants = ESConsts,
                        data = ESData, inits = ESInits)



###### testing ######



###### Compiling the model and MCMC ######

# nimbleMCMC compiles the model. If u want to customize your MCMC
# compile it yourself using:
#C.LipModel <- compileNimble(lipModel)

monitors <- c("b0","b1","b2")


cmodel  <- compileNimble(ESModel)
mcmc    <- buildMCMC(ESModel, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = ESModel)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)


predictive_check <- function(monitors, nSamp, n,
                             cmodel, samples, seed = 123,
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


y_posterior <- predictive_check(monitors = monitors,
                                nSamp = nrow(samples),
                                n = 78,
                                cmodel = cmodel,
                                samples = samples)

ES.sf$posteriorMean <- colMeans(y_posterior)

ggplot()+
  geom_sf(data = ES.sf, aes(fill = Y))

ggplot()+
  geom_sf(data = ES.sf, aes(fill = posteriorMean))

ES.sf$err <- ES.sf$Y - ES.sf$posteriorMean

ggplot()+
  geom_sf(data = ES.sf, aes(fill = err))

moran(ES.sf$err, lw, length(nb), Szero(lw))
moran.test(ES.sf$err,lw, alternative="two.sided")


###################

ESCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  b2 ~ dnorm(0,1)
  tau ~ dgamma(0.1, 0.1)
  
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau, zero_mean = 1)
  
  for (i in 1:N){
    
    log(theta[i]) <- b0 + b1*x1[i] + b2*x2[i] + s[i]
    y[i] ~ dpois(theta[i])
    
  }
})

ESAdj <- as.carAdjacency(adj)

ESConsts <- list(N = 78, L = length(ESAdj$adj),
                  adj = ESAdj$adj,
                  w = ESAdj$weights, num = ESAdj$num)
ESData <- list(y = ES.sf$Y, x1 = ES.sf$x1, x2 = ES.sf$x2)
ESInits <- list(b0 = 0, b1 = 0, b2 = 0, tau = 1)


ESModel <- nimbleModel(code = ESCode, name = "ES", constants = ESConsts,
                       data = ESData, inits = ESInits)



###### testing ######


###### Compiling the model and MCMC ######

# nimbleMCMC compiles the model. If u want to customize your MCMC
# compile it yourself using:
#C.LipModel <- compileNimble(lipModel)

monitors <- c("b0","b1","b2")


cmodel  <- compileNimble(ESModel)
mcmc    <- buildMCMC(ESModel, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = ESModel)
samples <- runMCMC(cmcmc, niter = 3000, nburnin = 1000)

y_posterior_CAR <- predictive_check(monitors = monitors,
                                nSamp = nrow(samples),
                                n = 78,
                                cmodel = cmodel,
                                samples = samples)


ES.sf$posteriorMean_CAR <- colMeans(y_posterior_CAR)

ES.sf$err_CAR <- ES.sf$Y - ES.sf$posteriorMean_CAR


ggplot()+
  geom_sf(data = ES.sf, aes(fill = posteriorMean_CAR))


ggplot()+
  geom_sf(data = ES.sf, aes(fill = err_CAR))

moran(ES.sf$err_CAR, lw, length(nb), Szero(lw))
moran.test(ES.sf$err_CAR,lw, alternative="two.sided")
