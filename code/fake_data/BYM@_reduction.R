library(ggplot2)
library(ggmap)
library(rgdal)
library(maps)
library(mapdata)
library(broom)
library(MASS)
library(sf)
library(sp)
library(spdep)
library(rwc)
library(dplyr)
library(nimble)
library(INLA)

set.seed(65465)

base_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\code\\maps"


ES.shp <- readOGR(
  dsn = paste0(base_dir, '\\shapefiles'),
  layer = 'ES_G'
)

ES.fortified <- tidy(ES.shp)
ES.sf <- st_as_sf(ES.shp)
ES.sf$G <- as.numeric(ES.sf$G)
N <- 78
nb <- poly2nb(ES.shp)
adj <- nb2mat(nb,   style = "B", zero.policy = TRUE)
lw <-  nb2listw(nb, style = "B", zero.policy=TRUE)

nb2INLA('nb', nb)
g <- inla.read.graph(filename = "nb")
#Create matrix Q
Q <- -inla.graph2matrix(g)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
# Scale intrinsic GMRF so that geometric mean of marginal variance is 1
Q.scaled <- inla.scale.model(Q,constr=list(A=matrix(1,1,N),e=0))
scale <- Q.scaled[1,1]/Q[1,1]

Preccc <- ginv(matrix(Q.scaled, ncol = 78, nrow = 78))

alec <- rmvnorm(n=1,mean=rep(0,78), sigma = ginv(matrix(Q.scaled, ncol = 78, nrow = 78)))

ndi <- rep(0,78)

for(i in 1:78){
  ndi[i] <- sum(adj[i,])
}

rho <- 0.99
Precisao <- diag(ndi) - rho*adj
CAR <- rnorm.Q(Precisao, zero.constraint = TRUE)
ES.sf$CAR <- as.vector(CAR)

ggplot()+
  geom_sf(data = ES.sf, aes(fill = G))+
  scale_fill_continuous(type="viridis", guide = 'none')+
  theme_void()


moran(ES.sf$G, lw, length(nb), Szero(lw))
moran.test(alec,lw, alternative="two.sided")

beta0 <- 0.4
beta1 <- 0.5
beta2 <- -0.5
beta3 <- 0.3

x1 <- rnorm(n=78, 0,1)
x2 <- rnorm(n=78, 0,1)

ES.sf$x1 <- x1
ES.sf$x2 <- x2

um.sobre.tau <- 0.5

log.mu <- beta0 + ES.sf$x1*beta1 + ES.sf$x2*beta2 + ES.sf$G*beta3 
log.mu_const <- beta0
log.mu_car <- beta0 + um.sobre.tau*alec
log.mu_orle <- beta0 + um.sobre.tau*rnorm(78,0,1)

mu <- exp(log.mu)
mu
y <- rpois(n=78,mu)

ES.sf$Y <- y

ES.sf$Y_const <- rpois(n=78,exp(log.mu_const))
ES.sf$Y_ORLE <- rpois(n=78,exp(log.mu_orle))
ES.sf$Y_CAR <- rpois(n=78,exp(log.mu_car))



BYM2.code <- nimbleCode({
  
  for (i in 1:N){
    # BYM2
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- b0 + bym[i]
    bym[i] <- (sqrt(1-rho)*u[i] + sqrt(rho/scale)*s[i])/sqrt(tau.b)
    u[i] ~ dnorm(0,1)
    
  }
  
  b0 ~ dnorm(0,1) #intercept prior
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], 1, zero_mean = 1) # CAR
  rho ~ dbeta(1,1) #Mixing parameter
  tau.b ~ dgamma(1, 0.1) #islands precision prior
  um.sobre.raiz.de.tau.b <- 1/sqrt(tau.b)
 
})

n.chains <- 1
n.burnin <- 10000
n.iter   <- 15000
ESAdj <- as.carAdjacency(adj)
ESConsts <- list(N = 78, L = length(ESAdj$adj),
                 adj = ESAdj$adj, scale = scale,
                 w = ESAdj$weights, num = ESAdj$num)

ESData_const <- list(y = ES.sf$Y_const)
ESData_CAR <- list(y = ES.sf$Y_CAR)
ESData_ORLE <- list(y = ES.sf$Y_ORLE)

ESInits <- list(b0 = 0, s = rnorm(78,0,1), rho = 0.5, tau.b = 1,
                u = rnorm(78,0,1))

monitors <- c("b0", "rho", "um.sobre.raiz.de.tau.b")

BYM2samples_const <- nimbleMCMC(code = BYM2.code,
                          data = ESData_const,
                          constants = ESConsts, 
                          inits = ESInits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          samplesAsCodaMCMC = TRUE, 
                          summary = TRUE, 
                          WAIC = TRUE)

BYM2samples_ORLE <- nimbleMCMC(code = BYM2.code,
                                data = ESData_ORLE,
                                constants = ESConsts, 
                                inits = ESInits,
                                monitors = monitors,
                                niter = n.iter,
                                nburnin = n.burnin,
                                nchains = n.chains, 
                                samplesAsCodaMCMC = TRUE, 
                                summary = TRUE, 
                                WAIC = TRUE)

BYM2samples_CAR <- nimbleMCMC(code = BYM2.code,
                                data = ESData_CAR,
                                constants = ESConsts, 
                                inits = ESInits,
                                monitors = monitors,
                                niter = n.iter,
                                nburnin = n.burnin,
                                nchains = n.chains, 
                                samplesAsCodaMCMC = TRUE, 
                                summary = TRUE, 
                                WAIC = TRUE)


BYM2samples_const$summary
BYM2samples_ORLE$summary
BYM2samples_CAR$summary
