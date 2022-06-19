library(nimble)
library(MASS)
library(sp)
library(SpatialEpi)
library(spdep)
library(ggplot2)


##################
####SIMULATION####
##################
simulation <- function(n.regions, n.coef, cov.mat,
                    mu.beta = 0, sd.beta = 1,
                    adj.mat, tau, g, sd,b0, boolCAR = FALSE){
  
  X <- mvrnorm(n.regions, mu = rep(0,n.coef), Sigma = cov.mat)
  beta <- rnorm(n.coef, mu.beta, sd.beta)
  err <- rnorm(n.regions, 0, sd)
  
  if(boolCAR){
    Adj <- as.carAdjacency(adj.mat)
    CAR <- rcar_proper(n=1,mu = rep(0,n.regions),adj =  Adj$adj,
                       num = Adj$num, tau = tau, gamma = g)
    
    log_theta <- b0 + X %*% beta + CAR + err
  }else{
    log_theta <- b0 + X %*% beta + err
  }
  
  
  Y <- rep(0,n.regions)
  
  for(i in 1:n.regions){
    Y[i] <- rpois(n=1, exp(log_theta[i]))
  }
  
  return(list(Y = Y, X = X, beta = beta))
  }

##################
####PARAMETERS####
##################

data(scotland)
d <- scotland$data
rownames(d) <- d$county
map <- scotland$spatial.polygon

map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)
nb <- poly2nb(map)
adj.mat <- nb2mat(nb, style = "B", zero.policy = TRUE)
adj.mat <- adj.mat[,-(which(colSums(adj.mat)==0))]
adj.mat <- adj.mat[ -(which(rowSums(adj.mat)==0)),]

Adj <- as.carAdjacency(adj.mat)

N <- length(adj.mat[,1])
p <- 10
cov.mat <- nimble::identityMatrix(p)
tau <- 1
g <- 0.01
sd <- 0.5
sd.beta <- 0.5
b0 <- 2



sim.noCAR <- simulation(n.regions = N, n.coef = p, cov.mat = identityMatrix(p),
        adj.mat = adj.mat, tau = tau, g = g, sd = sd, sd.beta = sd.beta,
        b0 = b0)

sim.CAR <- simulation(n.regions = N, n.coef = p, cov.mat = identityMatrix(p),
                     adj.mat = adj.mat, tau = tau, g = g, sd = sd, sd.beta = sd.beta,
                     b0 = b0, boolCAR = TRUE)
##################
######NIMBLE######
##################



CodeNoCAR <- nimbleCode({
  
  # priors
  sigma ~ dgamma(1, 1)
  for(i in 1:p+1){
    b[i] ~ dnorm(0,1)
  }
  
  for (i in 1:N){
    err[i] ~ dnorm(0,sigma)
    log(theta[i]) <- b[1] + inprod(X[i,], b[2:(p+1)]) + err[i]# + s[i]
    y[i] ~ dpois(theta[i])
  }
})

CodeCAR <- nimbleCode({
  
  # priors
  sigma ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  g ~ dgamma(1, 1)
  
  for(i in 1:p+1){
    b[i] ~ dnorm(0,1)
  }
  
  s[1:N] ~ dcar_proper(mu = mu[1:N], adj = adj[1:L], num = num[1:N],
  tau = tau, gamma = g)
  
  for (i in 1:N){
    mu[i] <- mu0
    err[i] ~ dnorm(0,sigma)
    log(theta[i]) <- b[1] + inprod(X[i,], b[2:(p+1)]) + err[i] + s[i]
    y[i] ~ dpois(theta[i])
  }
})

Consts <- list(N = N,p = p, L = length(Adj$adj),
                  adj = Adj$adj,
                  num = Adj$num,
               mu0=0)

Data <- list(X = sim.noCAR$X, y = sim.noCAR$Y)

Inits <- list(tau = 1, g = 0, sigma = 1, s = rep(0,N),b = rep(0,p+1))


##################
###sim. sem CAR###
##################

Model1 <- nimbleModel(code = CodeCAR, name = "fake", constants = Consts,
                        data = Data, inits = Inits)

Model2 <- nimbleModel(code = CodeNoCAR, name = "fake", constants = Consts,
                      data = Data, inits = Inits)

mcmc.out1 <- nimbleMCMC(code = CodeCAR, constants = Consts,
                       data = Data, inits = Inits,
                       nchains = 3, niter = 10000, summary = TRUE,
                       WAIC = TRUE, nburnin = 1000)

mcmc.out2 <- nimbleMCMC(code = CodeNoCAR, constants = Consts,
                        data = Data, inits = Inits,
                        nchains = 3, niter = 10000, summary = TRUE,
                        WAIC = TRUE, nburnin = 1000)


##################
###sim. com CAR###
##################

Data <- list(X = sim.CAR$X, y = sim.CAR$Y)

Model3 <- nimbleModel(code = CodeCAR, name = "fake", constants = Consts,
                      data = Data, inits = Inits)

Model4 <- nimbleModel(code = CodeNoCAR, name = "fake", constants = Consts,
                      data = Data, inits = Inits)

mcmc.out3 <- nimbleMCMC(code = CodeCAR, constants = Consts,
                        data = Data, inits = Inits,
                        nchains = 3, niter = 10000, summary = TRUE,
                        WAIC = TRUE, nburnin = 1000)

mcmc.out4 <- nimbleMCMC(code = CodeNoCAR, constants = Consts,
                        data = Data, inits = Inits,
                        nchains = 3, niter = 10000, summary = TRUE,
                        WAIC = TRUE, nburnin = 5000)


