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
  geom_sf(data = ebola_IN, aes(fill = log(counts)))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "logCasos")+
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



exclude <- c('location', 'OBJECTID', 'F1', 'location_1', 'order_',
             'country', 'include', 'counts', 'Offset_', 'Case_count',
             'Sequence_c', 'Pop_Size', 'SIR', 'area', 'perimetro',
             'geometry', 'locations')

colunas <- colnames(ebola_IN)
teste <- ebola_IN[colunas[!(colunas  %in% exclude)]]
teste
N.countries <- length(teste$geconMN)
N.covariates <- length(colnames(teste))-1
#nimble was not accepting ebola_IN.covariates, so had to find a quick fix
cov.ebola <- matrix(rep(0,N.countries*N.covariates),
                    nrow = N.countries, ncol = N.covariates)

for(i in 1:N.countries){
  for(j in 1:N.covariates){
    cov.ebola[i,j] <- as.numeric(teste[i,j])[1]
  }
}






#Simple model
modelo1_simples <- nimbleCode({
  
  for (i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p])
    
  }
  
  b0 ~ dnorm(0,0.001)
  beta[1:p] ~ dmnorm(mean = m[1:p], cov = cov[1:p, 1:p])
  
})


N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])

n.chains <- 2
n.burnin <- 200000
n.iter   <- 300000


monitors <- c("beta", "mu", "b0")

modelo1_simplesConsts <- list(N = N.countries,
                              p = N.covariates,
                              m = rep(0,N.covariates),
                              cov = identityMatrix(N.covariates))

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
modelo1_simplesInits <- list(b0 =  -1, beta = rnorm(N.covariates,0,1))


#Choose your way of running the model: compareMCMC or nimbleMCMC

comparr.model1 <- compareMCMCs( modelInfo = list(code = modelo1_simples,
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


waic.model1 <- nimbleMCMC(code = modelo1_simples,
                          constants = modelo1_simplesConsts,
                          data = modelo1_simplesData,
                          inits = modelo1_simplesInits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          thin = 10,
                          summary = TRUE,
                          WAIC = TRUE)


##################################################################################


#ORLE model
modelo2_ORLE <- nimbleCode({
  
  
  for (i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + u[i]
    u[i] ~ dnorm(0, tau.b)
    
  }
  
  tau.b ~ dgamma(1,0.01)
  b0 ~ dnorm(0,0.001) 
  beta[1:p] ~ dmnorm(mean = m[1:p], cov = cov[1:p, 1:p])
  
})


N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])

n.chains <- 2
n.burnin <- 200000
n.iter   <- 300000


monitors <- c("beta", "mu", "b0")

modelo2Consts <- list(N = N.countries,
                              p = N.covariates,
                              m = rep(0,N.covariates),
                              cov = identityMatrix(N.covariates)*100)


modelo2Data <- list(y = ebola_IN$counts, X = cov.ebola,
                            e = ebola_IN$Offset_/100000)

modelo2Inits <- list(b0 =  3, beta = rnorm(N.covariates,0,1),
                     u = rnorm(N.countries,0,1), tau.b = 1)

#Choose your way of running the model: compareMCMC or nimbleMCMC

comparr.model2 <- compareMCMCs( modelInfo = list(code = modelo2_ORLE,
                                             constants = modelo2Consts,
                                             data = modelo2Data,
                                             inits = modelo2Inits),
                            MCMCcontrol = list(niter = n.iter,
                                               thin = 10,
                                               burnin = n.burnin),
                            MCMCs = "nimble",
                            monitors = monitors,
                            metrics = c("mean", "median", "sd", "CI95_low",
                                        "CI95_upp", "efficiency_coda"))


waic.model2 <- nimbleMCMC(code = modelo2_ORLE,
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
#ORLE model
modelo3_BYM <- nimbleCode({
  
  
  for (i in 1:N){
    
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + u[i] + s[i]
    u[i] ~ dnorm(0, tau.b)
    
  }
  
  tau.s ~ dgamma(1,0.1)
  tau.b ~ dgamma(1,0.1)
  b0 ~ dnorm(0,0.001) 
  for(i in 1:p){
    beta[i] ~ dnorm(0,1)
  }
  #beta[1:p] ~ dmnorm(mean = m[1:p], cov = cov[1:p, 1:p])
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau.s, zero_mean = 1)
})


N.countries <- length(ebola_IN$country)
N.covariates <- length(ebola_IN.covariates[1,])

n.chains <- 2
n.burnin <- 200000
n.iter   <- 400000


monitors <- c("beta", "mu", "b0", "tau.s", "tau.b")

modelo3Consts <- list(N = N.countries,
                      p = N.covariates,
                      m = rep(0,N.covariates),
                      cov = identityMatrix(N.covariates),
                      L = length(EbolaAdj$adj),
                      adj = EbolaAdj$adj, w = EbolaAdj$weights, num = EbolaAdj$num)


modelo3Data <- list(y = ebola_IN$counts, X = cov.ebola,
                    e = ebola_IN$Offset_/100000)

modelo3Inits <- list(b0 =  0, beta = rep(0, N.covariates),
                     u = rnorm(N.countries,0,1), tau.b = 1,
                     s = rnorm(N.countries,0,1), tau.s = 1 )

#Choose your way of running the model: compareMCMC or nimbleMCMC

comparr.model3 <- compareMCMCs( modelInfo = list(code = modelo3_BYM,
                                                 constants = modelo3Consts,
                                                 data = modelo3Data,
                                                 inits = modelo3Inits),
                                MCMCcontrol = list(niter = n.iter,
                                                   thin = 20,
                                                   burnin = n.burnin),
                                MCMCs = "nimble",
                                monitors = monitors,
                                metrics = c("mean", "median", "sd", "CI95_low",
                                            "CI95_upp", "efficiency_coda"))


waic.model3 <- nimbleMCMC(code = modelo3_BYM,
                          constants = modelo3Consts,
                          data = modelo3Data,
                          inits = modelo3Inits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          thin = 10,
                          samples = TRUE,
                          summary = TRUE,
                          WAIC = TRUE)











######################
# Nimble model

BYM2Code <- nimbleCode({
  
  for (i in 1:N){
    # BYM2
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + inprod(X[i,1:p],beta[1:p]) + bym[i]
    bym[i] <- (sqrt(1-rho)*u[i] + sqrt(rho/scale)*s[i])/sqrt(tau.b)
    u[i] ~ dnorm(0,1)
    
  }
  
  b0 ~ dnorm(0,0.001) #intercept prior
  for(i in 1:p){
    beta[i] ~ dnorm(0, 1)
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
thinning <- 20


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
       u = rnorm(N.countries,0,1))

BYM2samples <- nimbleMCMC(code = BYM2Code,
                          data = BYM2Data,
                          constants = BYM2Consts, 
                          inits = BYM2Inits,
                          monitors = monitors,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          thin = thinning,
                          samplesAsCodaMCMC = TRUE, 
                          summary = TRUE, 
                          WAIC = TRUE)

#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(BYM2samples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.5) 

comparr.modelbym2 <- compareMCMCs( modelInfo = list(code = BYM2Code,
                                                 constants = BYM2Consts,
                                                 data = BYM2Data,
                                                 inits = BYM2Inits),
                                MCMCcontrol = list(niter = n.iter,
                                                   thin = thinning,
                                                   burnin = n.burnin),
                                MCMCs = "nimble",
                                monitors = monitors,
                                metrics = c("mean", "median", "sd", "CI95_low",
                                            "CI95_upp", "efficiency_coda"))

BYM2samples$summary$all.chains[c('gamma[1]'),]
BYM2samples$WAIC
