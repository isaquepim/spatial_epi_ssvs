library(nimble)
library(ggplot2)
library(MASS)
library(dplyr)
library(compareMCMCs)

set.seed(123)

#################################
######### SCENARIO I ############
#################################

n <- 200
p <- 7
beta <- c(3,0,0,-1,0,2,0)

err.sigma <- 3


err <- rnorm(n = n, sd = err.sigma)
Sigma <- matrix(rep(0.5,n*p),nrow=p,ncol=p)
diag(Sigma) <- 1



X <- mvrnorm(n,mu = rep(0,p), Sigma = Sigma)

y <- X %*% beta + err

#################################
######### THE MODEL #############
#################################


model.MC <- nimbleCode({
  
  for(i in 1:p){
    gamma[i] ~ dbern(0.5)
    mix[i] <- ((1-gamma[i])*tau + gamma[i]*c*tau)
    beta[i] ~ dnorm(mean = 0, sd = mix[i])
  }
  
  sigma2 ~ dinvgamma(shape = nu/2,scale = nu*lambda/2)
  
  for(i in 1:N){
    y[i] ~ dnorm(mean = inprod(X[i,1:p],beta[1:p]),sd = sqrt(sigma2))
  }
  
})


tau <- 0.01
c <- 5000

nu <- 10
lambda <- 1

model.consts <- list(N = n, p = p, tau = tau, c = c,
                     nu = nu, lambda = lambda)

model.data <- list(X = X, y = c(y))
model.inits <- list(beta = rnorm(p,0,1), gamma = rep(1,p),
                    sigma2 = 1)
model.Model <- nimbleModel(code = model.MC,
                           name = "Regression case I",
                           constants = model.consts,
                           data = model.data,
                           inits = model.inits)

monitors <- c("beta", "gamma", "sigma2")

cmodel  <- compileNimble(model.Model)
mcmc    <- buildMCMC(model.Model, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = model.Model)
samples_MC <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)

samplesSummary(samples_MC)

df <-  as.data.frame(samples_MC[,c('beta[3]', 'gamma[3]')])
colnames(df) <- c('x','g')

ggplot()+
  geom_density(data = subset(df, g == 0), aes(x = x), fill = 'red', alpha = 0.3)


model.KM <- nimbleCode({
  
  for(i in 1:p){
    gamma[i] ~ dbern(0.5)
    b[i] ~ dnorm(0,1)
    beta[i] <- gamma[i]*b[i]
  }
  
  sigma2 ~ dinvgamma(shape = nu/2,scale = nu*lambda/2)
  
  for(i in 1:N){
    y[i] ~ dnorm(mean = inprod(X[i,1:p],beta[1:p]),sd = sqrt(sigma2))
  }
  
})

nu <- 10
lambda <- 1

model.consts <- list(N = n, p = p,
                     nu = nu, lambda = lambda)

model.data <- list(X = X, y = c(y))
model.inits <- list(b = rnorm(p,0,1), gamma = rep(1,p),
                    sigma2 = 1)
model.Model <- nimbleModel(code = model.KM,
                           name = "Regression case I",
                           constants = model.consts,
                           data = model.data,
                           inits = model.inits)

monitors <- c("beta", "gamma", "sigma2")

cmodel  <- compileNimble(model.Model)
mcmc    <- buildMCMC(model.Model, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = model.Model)
samples_KM <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)

samplesSummary(samples_KM)

df <-  as.data.frame(samples_KM[,c('beta[4]', 'gamma[4]')])
colnames(df) <- c('x','g')

ggplot()+
  geom_density(data = df, aes(x = x), fill = 'red', alpha = 0.3)


comparr <- compareMCMCs( modelInfo = list(code = model.KM,
                                          constants = model.consts,
                                          data = model.data,
                                          inits = model.inits),
                            MCMCcontrol = list(niter = 10000, burnin = 2000),
                            MCMCs = "nimble",
                            monitors = monitors,
                            metrics = c("mean", "median", "sd", "CI95_low",
                                        "CI95_upp", "efficiency_coda"))
comparr.MC <- compareMCMCs( modelInfo = list(code = model.MC,
                                             constants = model.consts,
                                             data = model.data,
                                             inits = model.inits),
                            MCMCcontrol = list(niter = 10000, burnin = 2000),
                            MCMCs = "nimble",
                            monitors = monitors,
                            metrics = c("mean", "median", "sd", "CI95_low",
                                        "CI95_upp", "efficiency_coda"))

make_MCMC_comparison_pages(comparr, dir = getwd())
