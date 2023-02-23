library(nimble)
library(spdep)
library(sp)
library(ggplot2)
library(INLA)
library(coda)
library(sf)
library(compareMCMCs)


##############################################
#LOADING SHAPEFILE

ebolain_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\data\\ebola_IN.shp"

ebola_IN <- read_sf(ebolain_dir)

ebola_IN$SIR <- ebola_IN$counts/ebola_IN$Pop_Size
ebola_IN$PopDens <- ebola_IN$Pop_Size/ebola_IN$area

nb.ebola <- poly2nb(ebola_IN)
W <- nb2mat(nb.ebola, style = "B", zero.policy = TRUE)
EbolaAdj <- as.carAdjacency(W)

ggplot()+
  geom_sf(data = ebola, aes(fill = include))+
  scale_fill_manual(values=c(OUT = "#74736E", IN = "#FAAC5A"))+
  labs(fill="Incluir")+
  theme_void()

ggplot()+
  geom_sf(data = ebola_IN, aes(fill = log(pm)))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "Casos")+
  labs(title = "ORLE")+
  theme_void()

ggplot()+
  geom_sf(data = ebola_IN, aes(fill = log(counts)))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "Casos")+
  labs(title = "DADO")+
  theme_void()
ggplot()+
  geom_sf(data = ebola_IN, aes(fill = Offset_))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "Casos")+
  theme_void()

ggplot()+
  geom_sf(data = ebola_IN, aes(fill = log(bym2pm)))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "Casos", )+
  labs(title = "BYM2")+
  theme_void()

ebola_IN$pm <- posterior_means$median

#####################
#neighbour list

colnames(ebola_IN)

######################
# subsetting variables

variables.to.include <- c("location","country","include", "counts", "Offset_",
                          "Case_count" ,"Sequence_c" ,"Pop_Size",  "geconMN",
                          "geconMIN", "geconMAX","geconSTD","pdensMN",        
                          "tt50kMN" , "tt100kMN"      , "tt500kMN", "altMN" , "tempMN"  ,
                          "tmpssMN" , "precMN", "precssMN"  , "Introadmin", "Introdista")

covariates <- c("geconMN", "geconMIN", "geconMAX","geconSTD","pdensMN",        
                "tt50kMN" , "tt100kMN"  , "tt500kMN", "altMN" , "tempMN"  ,
                "tmpssMN" , "precMN", "precssMN"  , "Introadmin", "Introdista")

ebola_IN.filtered <- ebola_IN[variables.to.include]
ebola_IN.matrix <- as.matrix(ebola_IN.filtered)
ebola_IN.covariates <- ebola_IN.matrix[,covariates]
ebola_IN.covariates




#Simple model
modelo1_simples_ssvs <- nimbleCode({
  
  for (i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p])
    
  }
  
  b0 ~ dnorm(0,0.01)
  for(i in 1:p){
    b[i] ~ dnorm(0,0.1)
    gamma[i] ~ dbern(q)
    beta[i] <- b[i]*gamma[i]
  }
})

N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])

n.burnin <- 200000
n.iter   <- 300000


monitors <- c("beta", "mu", "b0")

modelo1_simplesConsts <- list(N = N.countries,
                              p = N.covariates,
                              m = rep(0,N.covariates),
                              cov = identityMatrix(N.covariates),
                              q = 0.5)

cov.ebola <- matrix(rep(0,N.countries*N.covariates),
                    nrow = N.countries, ncol = N.covariates)

for(i in 1:N.countries){
  for(j in 1:N.covariates){
    cov.ebola[i,j] <- as.numeric(ebola_IN.covariates[i,j])
  }
}

modelo1_simplesData <- list(y = ebola_IN$counts, X = cov.ebola,
                            e = ebola_IN$Offset_/100000)

# initiate chain from 2 differente places
modelo1_simplesInits <- list(b0 =  -1, b = rnorm(N.covariates,0,1),
                             gamma = rep(1,N.covariates))


#Choose your way of running the model: compareMCMC or nimbleMCMC

comparr.model1_ssvs <- compareMCMCs( modelInfo = list(code = modelo1_simples_ssvs,
                                                 constants = modelo1_simplesConsts,
                                                 data = modelo1_simplesData,
                                                 inits = modelo1_simplesInits),
                                MCMCcontrol = list(niter = n.iter,
                                                   thin = 10,
                                                   burnin = n.burnin),
                                MCMCs = "nimble",
                                monitors = monitors,
                                metrics = c("mean", "median", "sd", "CI95_low",
                                            "CI95_upp", "efficiency_coda"))


waic.model1_ssvs <- nimbleMCMC(code = modelo1_simples_ssvs,
                          constants = modelo1_simplesConsts,
                          data = modelo1_simplesData,
                          inits = modelo1_simplesInits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          thin = 10,
                          summary = TRUE,
                          WAIC = TRUE)


#RMSE
sqrt(mean((colMeans(waic.model1_ssvs$samples[,mapply(function(x) paste0('mu[',x,']'), c(1:63))]) - ebola_IN$counts)^2))


##################################################################################


#ORLE model
modelo2_ORLE_ssvs <- nimbleCode({
  
  
  for (i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + u[i]
    u[i] ~ dnorm(0, tau.b)
    
  }
  
  tau.b ~ dgamma(1,0.01)
  b0 ~ dnorm(0,0.001) 
  for(i in 1:p){
    b[i] ~ dnorm(0,0.1)
    gamma[i] ~ dbern(q)
    beta[i] <- b[i]*gamma[i]
  }
  
})


N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])

n.chains <- 2
n.burnin <- 200000
n.iter   <- 400000


monitors <- c("beta", "mu", "b0")

modelo2Consts <- list(N = N.countries,
                      p = N.covariates,
                      m = rep(0,N.covariates),
                      cov = identityMatrix(N.covariates)*100,
                      q = 0.5)


modelo2Data <- list(y = ebola_IN$counts, X = cov.ebola,
                    e = ebola_IN$Offset_/100000)

modelo2Inits <- list(b0 =  3, beta = rnorm(N.covariates,0,1),
                     u = rnorm(N.countries,0,1), tau.b = 1,
                     gamma = rep(1,N.covariates))

#Choose your way of running the model: compareMCMC or nimbleMCMC

comparr.model2_ssvs <- compareMCMCs( modelInfo = list(code = modelo2_ORLE_ssvs,
                                                 constants = modelo2Consts,
                                                 data = modelo2Data,
                                                 inits = modelo2Inits),
                                MCMCcontrol = list(niter = n.iter,
                                                   thin = 20,
                                                   burnin = n.burnin),
                                MCMCs = "nimble",
                                monitors = monitors,
                                metrics = c("mean", "median", "sd", "CI95_low",
                                            "CI95_upp", "efficiency_coda"))


waic.model2_ssvs <- nimbleMCMC(code = modelo2_ORLE_ssvs,
                          constants = modelo2Consts,
                          data = modelo2Data,
                          inits = modelo2Inits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          thin = 10,
                          samples = TRUE,
                          summary = TRUE,
                          WAIC = TRUE)
##################################################################################

posterior_means[mapply(function(x) paste0('mu[',x,']'), c(1:63)),]



######################
# Nimble model

BYM2Code_ssvs <- nimbleCode({
  
  for (i in 1:N){
    # BYM2
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + bym[i]
    bym[i] <- (sqrt(1-rho)*u[i] + sqrt(rho/scale)*s[i])/sqrt(tau.b)
    u[i] ~ dnorm(0,1)
    
  }
  
  b0 ~ dnorm(0,0.001) #intercept prior
  for(i in 1:p){
    b[i] ~ dnorm(0,0.1)
    gamma[i] ~ dbern(q)
    beta[i] <- b[i]*gamma[i]
  }
  
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], 1, zero_mean = 1) # CAR
  rho ~ dbeta(1,1) #Mixing parameter
  
  tau.b ~ dgamma(1, 0.1) #islands precision prior
  
})


N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])


nb.ebola <- poly2nb(ebola_IN)
W <- nb2mat(nb.ebola, style = "B", zero.policy = TRUE)
EbolaAdj <- as.carAdjacency(W)

nb2INLA('nb_ebola', nb.ebola)
g <- inla.read.graph(filename = "nb_ebola")
#Create matrix Q
Q <- -inla.graph2matrix(g)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
# Scale intrinsic GMRF so that geometric mean of marginal variance is 1
Q.scaled <- inla.scale.model(Q,constr=list(A=matrix(1,1,N.countries),e=0))
scale <- Q.scaled[1,1]/Q[1,1]


EbolaAdj <- as.carAdjacency(W)

#mcmc

n.chains <- 1
n.burnin <- 200000
n.iter   <- 400000


monitors <- c("beta", "rho", "mu")

BYM2Consts <- list(N = N.countries,
                   p = N.covariates,
                   L = length(EbolaAdj$adj),
                   adj = EbolaAdj$adj, w = EbolaAdj$weights, num = EbolaAdj$num,
                   scale = scale)


BYM2Data <- list(y = ebola_IN$counts, X = cov.ebola,
                 e = ebola_IN$Offset_/100000)

# initiate chain from 4 differente places
BYM2Inits <- list(b0 =  0, beta = rnorm(N.covariates,0,1), rho = 0.5,
                  s = rnorm(N.countries,0,1), tau.b = 1,
                  u = rnorm(N.countries,0,1), gamma = rep(1,N.covariates))

BYM2samples_ssvs <- nimbleMCMC(code = BYM2Code_ssvs,
                          data = BYM2Data,
                          constants = BYM2Consts, 
                          inits = BYM2Inits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          thin = 20,
                          samplesAsCodaMCMC = TRUE, 
                          summary = TRUE, 
                          WAIC = TRUE)

#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(BYM2samples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.5) 



BYM2samples$summary$all.chains[c('gamma[1]'),]
BYM2samples$WAIC
