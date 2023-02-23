library(nimble)
library(coda)
library(sf)



ebola_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\data\\ebola_IN.shp"
ebola <- read_sf(ebola_dir)

covariates <- c("geconMN", "geconMIN", "geconMAX","geconSTD","pdensMN",        
                "tt50kMN" , "tt100kMN"  , "tt500kMN", "altMN" , "tempMN"  , "tmpssMN" , "precMN", "precssMN"  , "Introadmin", "Introdista")


N.countries <- length(ebola$country)
N.covariates <- length(covariates)

ebola.matrix <- as.matrix(ebola)
ebola.covariates <- ebola.matrix[,covariates]
ebola.covariates <- matrix(as.numeric(ebola.covariates),
                           nrow = N.countries, ncol = N.covariates)




ORLE.code <- nimbleCode({
  
  for(i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + ORLE[i]
    ORLE[i] ~ dnorm(0,tau.u)
    
  }
  
  b0 ~ dnorm(0, 0.001)
  for(i in 1:p){
    beta[i] ~ dnorm(0,0.001)
  }
  
  tau.u ~ dgamma(0.001,0.001)
  
})


#mcmc

n.chains <- 10
n.burnin <- 2000
n.iter   <- 100000
thinning <- 1

monitors <- c("b0", "beta", "tau.u")


ORLE.consts <- list(p = N.covariates, N = N.countries)
ORLE.data <- list(X = ebola.covariates, y = ebola$counts,
                  e = ebola$Offset_/100000)

ORLE.inits <- list()
for(i in 1:n.chains){
  new_init <- list(b0 = rnorm(1,0,1), beta = rnorm(N.covariates,0,1),
                   ORLE = rnorm(N.countries,0,1), tau.u = abs(rnorm(1,0,4)))
  ORLE.inits <- append(ORLE.inits, new_init)
}




ORLE.samples <- nimbleMCMC(code = ORLE.code,
                          data = ORLE.data,
                          constants = ORLE.consts, 
                          inits = ORLE.inits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains, 
                          thin = thinning,
                          samplesAsCodaMCMC = TRUE, 
                          summary = TRUE, 
                          WAIC = TRUE)


#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(ORLE.samples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.05) 

GR.diag$psrf[,"Point est."]





