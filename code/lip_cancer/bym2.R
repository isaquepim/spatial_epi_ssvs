library(nimble)
library(spdep)
library(SpatialEpi)
library(sp)
library(ggplot2)
library(INLA)


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

nb.no_islands <- subset(nb, card(nb) != 0)
W.no_islands <- W[-islands, -islands]
cases.no_islands <- d$cases[-islands]
aff.no_islands <- d$AFF[-islands]
E.no_islands <- d$expected[-islands]

cases.islands <- d$cases[islands]
aff.islands <- d$AFF[islands]
E.islands <- d$expected[islands]
#Calculating scale
N.no_islands <- length(cases.no_islands)
N.islands <- length(islands)

nb2INLA('nb', nb.no_islands)
g <- inla.read.graph(filename = "nb")
#Create matrix Q
Q <- -inla.graph2matrix(g)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
# Scale intrinsic GMRF so that geometric mean of marginal variance is 1
Q.scaled <- inla.scale.model(Q,constr=list(A=matrix(1,1,N.no_islands),e=0))
scale <- Q.scaled[1,1]/Q[1,1]
scale2 <- exp((1/nrow(Q))*sum(log(1/diag(Q.scaled))))
scale


#matrix to adjacency list
LipAdj <- as.carAdjacency(W.no_islands)



lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
  tau.err ~ dgamma(3.2,1.8)
  tau ~ dgamma(1, 1)
  rho ~ dbeta(1,2)
  
  
  s[1:Ni] ~ dcar_normal(adj[1:L], w[1:L], num[1:Ni], 1, zero_mean = 1)
  
  for (i in 1:Ni){
    
    err.i[i] ~ dnorm(0,tau.err)
    bym[i] <- (sqrt(1-rho)*err.i[i] + sqrt(rho/scale)*s[i])/sqrt(tau)
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
                  I = N.islands, scale = scale2)

lipData <- list(y = cases.islands, aff = aff.islands, e = E.islands,
                y.i = cases.no_islands, aff.i = aff.no_islands, e.i = E.no_islands)

lipInits <- list(b0 = 0, b1 = 0, tau = 1, tau.err = 1, s = rep(0,N.no_islands),
                 err = rnorm(E.islands, 0,1), err.i = rnorm(N.no_islands,0,1),
                 rho = 0.5)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)


###### Compiling the model and MCMC ######
monitors <- c("b0","b1","tau","tau.err", "rho")

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

library("devtools")
devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE)
