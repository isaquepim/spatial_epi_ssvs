# Reproducing George & McCulloch
# Canonical regression setup, large p
# Setup IV
# tuple = (10,500), indiference prior.


library(nimble)
library(ggplot2)
library(MASS)
library(dplyr)

set.seed(123)

#################################
######### SCENARIO I ############
#################################

n <- 120
p <- 30
beta <- c(rep(3,15),rep(0,15))

err.sigma <- 25


err <- rnorm(n = n, sd = err.sigma)

X_ <- mvrnorm(n,mu = rep(0,p), Sigma = identityMatrix(p))
Z <- rnorm(n = n, mean = 0, sd = 1)

X <- X_ + matrix(rep(Z,p), nrow = n, ncol = p)
y <- X %*% beta + err

#################################
######### THE MODEL #############
#################################


model.code <- nimbleCode({
  
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
model.inits <- list(beta = rep(0,p), gamma = rep(1,p),
                    sigma2 = 1)



model.Model <- nimbleModel(code = model.code,
                           name = "Regression case I",
                           constants = model.consts,
                           data = model.data,
                           inits = model.inits)

monitors <- c("beta", "gamma", "sigma2")

cmodel  <- compileNimble(model.Model)
mcmc    <- buildMCMC(model.Model, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = model.Model)
samples_3 <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)

samplesSummary(samples_3)

df <-  as.data.frame(samples_3[,c('beta[26]', 'gamma[26]')])
colnames(df) <- c('x','g')

ggplot()+
  geom_density(data = df, aes(x = x), fill = 'red', alpha = 0.3)

