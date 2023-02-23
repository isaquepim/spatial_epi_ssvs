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
p <- 60
beta <- c(rep(0,15),rep(1,15),rep(2,15),rep(-1,15))

err.sigma <- 1


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
    mix[i] <- ((1-gamma[i])*tau + gamma[i]*c*tau)
    beta[i] ~ dnorm(mean = 0, sd = mix[i])
  }
  
  
  sigma2 ~ dinvgamma(shape = nu/2,scale = nu*lambda/2)
  
  for(i in 1:N){
    y[i] ~ dnorm(mean = inprod(X[i,1:p],beta[1:p]),sd = sqrt(sigma2))
  }
})

# OLS estimate:

tau <- 0.01
c <- 500

nu <- 10
lambda <- 1

model.consts <- list(N = n, p = p, tau = tau, c = c,
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
samples_4 <- runMCMC(cmcmc, niter = 10000, nburnin = 2000)

samplesSummary(samples_4)

df <-  as.data.frame(samples_4[,c('beta[23]', 'gamma[23]')])
colnames(df) <- c('x','g')

ggplot()+
  geom_density(data = df, aes(x = x), fill = 'red', alpha = 0.3)





ggplot()+
  geom_line(data = subset(df, g == 0), aes(x = c(1:length(x)),y = x))

ggplot()+
  geom_line(data = subset(df, g == 1), aes(x = c(1:length(x)),y = x))

