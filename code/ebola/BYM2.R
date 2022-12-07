library(nimble)
library(coda)
library(INLA)
library(sf)
library(spdep)
library(corrplot)


ebola_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\data\\ebola_IN.shp"
ebola <- read_sf(ebola_dir)

covariates <- c("geconMN", "geconMIN", "geconMAX","geconSTD","pdensMN",        
                "tt50kMN" , "tt100kMN"  , "tt500kMN", "altMN" , "tempMN"  , "tmpssMN" , "precMN", "precssMN"  , "Introadmin", "Introdista")

covariates.nocor <- covariates[-c(1,2,7,8,9,11)]

N.countries <- length(ebola$country)
N.covariates <- length(covariates)

ebola.matrix <- as.matrix(ebola)
ebola.covariates <- ebola.matrix[,covariates]
ebola.covariates <- matrix(as.numeric(ebola.covariates),
                           nrow = N.countries, ncol = N.covariates)


corrplot(cor(ebola.covariates),method = 'number')





BYM2.code <- nimbleCode({
  
  for(i in 1:N){
    
    y[i] ~ dpois(mu[i])
    
    log(mu[i]) <- log(e[i]) + b0 + sum(X[i,1:p]*beta[1:p]) + bym[i]
    
    bym[i] <- (sqrt(1-rho)*u[i] + sqrt(rho/scale)*s[i])/sqrt(tau.r)
    u[i] ~ dnorm(0,1)
    
  }
  
  b0 ~ dnorm(0, 0.01)
  for(i in 1:p){
    beta[i] ~ dnorm(0,0.01)
  }
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], 1, zero_mean = 1) #CAR
  rho ~ dbeta(1,1) #Mixing parameter
  
  tau.r ~ dgamma(1,0.01)
  
})

nb.ebola <- poly2nb(ebola)
W <- nb2mat(nb.ebola, style = "B", zero.policy = TRUE)

nb2INLA('nb_ebola', nb.ebola)
g <- inla.read.graph(filename = "nb_ebola")
#Create matrix Q
Q <- -inla.graph2matrix(g)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
# Scale intrinsic GMRF so that geometric mean of marginal variance is 1
Q.scaled <- inla.scale.model(Q,constr=list(A=matrix(1,1,N.countries),e=0))
scale <- exp((1/nrow(Q))*sum(log(1/diag(Q.scaled))))

EbolaAdj <- as.carAdjacency(W)


#mcmc

n.chains <- 10
n.burnin <- 30000
n.iter   <- 6000000
thinning <- 1000

monitors <- c("b0", "beta", "tau.r", "rho")


BYM2.consts <- list(p = N.covariates, N = N.countries,
                    L = length(EbolaAdj$adj),
                    adj = EbolaAdj$adj, w = EbolaAdj$weights, num = EbolaAdj$num,
                    scale = scale)
BYM2.data <- list(X = ebola.covariates, y = ebola$counts,
                  e = ebola$Offset_/10000)

BYM2.inits <- list()
for(i in 1:n.chains){
  new_init <- list(b0 = rnorm(1,0,2), beta = rnorm(N.covariates,0,1),
                   rho = rbeta(1,1,1), 
                   u = rnorm(N.countries,0,2),
                   s = rnorm(N.countries,0,2),
                   tau.r = abs(rnorm(1,0,4)))
  BYM2.inits <- append(BYM2.inits, new_init)
}




BYM2.samples <- nimbleMCMC(code = BYM2.code,
                           data = BYM2.data,
                           constants = BYM2.consts, 
                           inits = BYM2.inits,
                           monitors = monitors,
                           niter = n.iter,
                           nburnin = n.burnin,
                           nchains = n.chains, 
                           thin = thinning,
                           samplesAsCodaMCMC = TRUE, 
                           summary = TRUE, 
                           WAIC = TRUE)


#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(BYM2.samples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.05) 

GR.diag$psrf[,"Point est."]


library(ggplot2)

p <- ggplot()

for(i in 1:n.chains){
  p <- p +  geom_line(data = as.data.frame(as.matrix(BYM2.samples$samples[paste0("chain",i)]))
                 , aes(x = 1:5970, y = `beta[15]`))
}

p


