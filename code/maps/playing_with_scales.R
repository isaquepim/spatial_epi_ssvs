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
library(INLA)
library(rwc)
library(dplyr)

library(nimble)

set.seed(65465)

base_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\code\\maps"


ES.shp <- readOGR(
  dsn = paste0(base_dir, '\\shapefiles'),
  layer = 'ES_G_regions'
)

ES.fortified <- tidy(ES.shp)
ES.sf <- st_as_sf(ES.shp)


nb <- poly2nb(ES.shp)
adj <- nb2mat(nb,   style = "B", zero.policy = TRUE)
lw <-  nb2listw(nb, style = "B", zero.policy=TRUE)

ES.sf$G <- as.numeric(ES.sf$G)

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


###################
#Dean



ESCode <- nimbleCode({
  for (i in 1:N){
    y[i] ~ dpois(theta[i])
    log(theta[i]) <- b0 + b1*x1[i] + b2*x2[i] + s[i]
  }
  # priors
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  b2 ~ dnorm(0,1)
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau, zero_mean = 1)
  
})

ES.sf1 <- subset(ES.sf, regiao == 1)

nb <- poly2nb(ES.sf)
adj <- nb2mat(nb,   style = "B", zero.policy = TRUE)
lw <-  nb2listw(nb, style = "B", zero.policy=TRUE)

W <- nb2mat(nb, style = "B", zero.policy = TRUE)

nb2INLA('nb', nb)
g <- inla.read.graph(filename = "nb")
#Create matrix Q
Q <- -inla.graph2matrix(g)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
# Scale intrinsic GMRF so that geometric mean of marginal variance is 1
Q.scaled <- inla.scale.model(Q,constr=list(A=matrix(1,1,78),e=0))
scale <- Q.scaled[1,1]/Q[1,1]

ESAdj <- as.carAdjacency(adj)

ESConsts <- list(N = 78, L = length(ESAdj$adj),
                 adj = ESAdj$adj,
                 tau = 0.001,
                 scale=scale2,
                 w = ESAdj$weights, num = ESAdj$num)

ESData <- list(y = ES.sf$Y, x1 = ES.sf$x1, x2 = ES.sf$x2)

ESInits <- list(b0 = 0, b1 = 0, b2 = 0,
                s = rnorm(78,0,1))


ESModel <- nimbleModel(code = ESCode, name = "ES", constants = ESConsts,
                       data = ESData, inits = ESInits)


monitors <- c("s")


cmodel  <- compileNimble(ESModel)
mcmc    <- buildMCMC(ESModel, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = ESModel)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 1000)


empirical_sd <- rep(0,78)

for(i in 1:78){
  empirical_sd[i] <- var(samples[,i])
}

sigma_gv <- exp(sum(log(empirical_sd))/78)
sigma_gv
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
