library(nimble)
library(ggplot2)
library(MASS)

set.seed(123)

###################
#######CASO1#######
#######||||########
#######VVVV########

#Caso 1 - Sem correlação, n > p

n1 <- 500
p1 <- 2

Sigma1 <- nimble::identityMatrix(p1)
err.sd1 <- 0.1

X1 <- mvrnorm(n = n1,mu = rep(0,p1) ,Sigma = Sigma1)

beta1 <- c(2,0)
#beta1 <- rnorm(n = p1, mean = 0, sd = 1)
#beta1[11:20] <- 0

b0 <- 1

err <- rnorm(n = n1, mean = 0, sd = err.sd1)

theta1 <- X1 %*% beta1 + err

y1 <- rep(0,n1)

for(i in 1:n1){
  y1[i] <- rpois(n = 1, exp(theta1[i] + b0))
}

# Nimble model

caso1 <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  
  for(i in 1:P){
    b[i] ~ dnorm(0,1)
    z[i] ~ dbern(0.5)
    bz[i] <- b[i]*z[i]
  }
  
  s ~ dinvgamma(10,1)
  
  for(i in 1:N){
    
    err[i] ~ dnorm(0,s)
    log(theta[i]) <- b0 + inprod(X[i,], bz[1:P]) + err[i]
    y[i] ~ dpois(theta[i])
  }
})


Const1 <- list(N = n1, P = p1)

Data1 <- list(y = y1, X = X1)

Inits1 <- list(z = rep(1,p1), b = rep(0,p1), s = 1, b0 = 0)


Model1 <- nimbleModel(code = caso1, name = "sns",
                        constants = Const1,
                        data = Data1, inits = Inits1)


monitors <- c("b", "z", "s", "b0")


cmodel  <- compileNimble(Model1)
mcmc    <- buildMCMC(Model1, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = Model1)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 5000)


samplesSummary(samples)


###################
#######CASO2#######
#######||||########
#######VVVV########

n2 <- 500
p2 <- 10

Sigma2 <- nimble::identityMatrix(p2)
err.sd2 <- 0.1

X2 <- mvrnorm(n = n2,mu = rep(0,p2) ,Sigma = Sigma2)

beta2 <- c(3,0,3,0,3,0,3,0,3,0)

b0 <- 1

err2 <- rnorm(n = n2, mean = 0, sd = err.sd2)

theta2 <- X2 %*% beta2 + err2

y2 <- rep(0,n2)

for(i in 1:n2){
  y2[i] <- rpois(n = 1, exp(theta2[i] + b0))
}

# Nimble model

caso2 <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  
  for(i in 1:P){
    b[i] ~ dnorm(0,1)
    z[i] ~ dbern(0.5)
    bz[i] <- b[i]*z[i]
  }
  
  s ~ dinvgamma(10,1)
  
  for(i in 1:N){
    
    err[i] ~ dnorm(0,s)
    log(theta[i]) <- b0 + inprod(X[i,], bz[1:P]) + err[i]
    y[i] ~ dpois(theta[i])
  }
})


Const2 <- list(N = n2, P = p2)

Data2 <- list(y = y2, X = X2)

Inits2 <- list(z = rep(1,p2), b = rep(0,p2), s = 1, b0 = 0)


Model2 <- nimbleModel(code = caso2, name = "sns",
                      constants = Const2,
                      data = Data2, inits = Inits2)


monitors <- c("b", "z", "s", "b0")


cmodel  <- compileNimble(Model2)
mcmc    <- buildMCMC(Model2, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = Model2)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 5000)


samplesSummary(samples)

###################
#######CASO3#######
#######||||########
#######VVVV########


#Normalizando as taxas para o intervalo (-0.7, 4.0)


n3 <- 500
p3 <- 10

Sigma3 <- nimble::identityMatrix(p3)
err.sd3 <- 0.1

X3 <- mvrnorm(n = n3,mu = rep(0,p3) ,Sigma = Sigma3)

beta3 <- c(3,0,3,0,3,0,3,0,3,0)

b0 <- 1

err3 <- rnorm(n = n3, mean = 0, sd = err.sd3)

theta3 <- b0 + X3 %*% beta3 + err3

theta3 <- (theta3 - mean(theta3))/sd(theta3)*0.875 + 2.35


y3 <- rep(0,n3)

for(i in 1:n3){
  y3[i] <- rpois(n = 1, exp(theta3[i]))
}

# Nimble model

caso3 <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  
  for(i in 1:P){
    b[i] ~ dnorm(0,1)
    z[i] ~ dbern(0.5)
    bz[i] <- b[i]*z[i]
  }
  
  s ~ dinvgamma(10,1)
  
  for(i in 1:N){
    
    err[i] ~ dnorm(0,s)
    log(theta[i]) <- b0 + inprod(X[i,], bz[1:P]) + err[i]
    y[i] ~ dpois(theta[i])
  }
})


Const3 <- list(N = n3, P = p3)

Data3 <- list(y = y3, X = X3)

Inits3 <- list(z = rep(1,p3), b = rep(0,p3), s = 1, b0 = 0)


Model3 <- nimbleModel(code = caso3, name = "sns",
                      constants = Const3,
                      data = Data3, inits = Inits3)


monitors <- c("b", "z", "s", "b0")


cmodel  <- compileNimble(Model3)
mcmc    <- buildMCMC(Model3, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = Model3)
samples <- runMCMC(cmcmc, niter = 10000, nburnin = 5000)


samplesSummary(samples)

















