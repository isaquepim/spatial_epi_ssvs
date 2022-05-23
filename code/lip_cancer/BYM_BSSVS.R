library(nimble)
library(sp)
library(SpatialEpi)
library(spdep)
library(ggplot2)

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
LipAdj <- nb2mat(nb, style = "B", zero.policy = TRUE)
non.zero <- LipAdj[,-(which(colSums(LipAdj)==0))]
non.zero <- non.zero[-(which(rowSums(non.zero)==0)),]
non.zero.adj <- as.carAdjacency(non.zero)
#matrix to adjacency list
LipAdj <- as.carAdjacency(LipAdj)




lipCode <- nimbleCode({
  
  # priors
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(0,1)
  sigma ~ dgamma(1, 1)
  tau ~ dgamma(1, 1)
  
  psi ~ dunif(0,1)
  z0 ~ dbern(psi)
  z1 ~ dbern(psi)
  
  zb0 <- b0*z0
  zb1 <- b1*z1
  
  # stand-alone priors for comparison
  p_b0 ~ dnorm(0,1)
  p_b1 ~ dnorm(0,1)
  p_sigma ~ dgamma(1, 1)
  p_tau ~ dgamma(1, 1)
  p_psi ~ dunif(0,1)
  
  
  s[1:N] ~ dcar_normal(adj[1:L], w[1:L], num[1:N], tau, zero_mean = 1)
  
  for (i in 1:N){
    
    err[i] ~ dnorm(0,sigma)
    log(theta[i]) <- zb0 + zb1*aff[i] + s[i] + err[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
})


###### Model ######
lipConsts <- list(N = length(d$cases), L = length(LipAdj$adj),
                  adj = LipAdj$adj,
                  w = LipAdj$weights, num = LipAdj$num)

lipData <- list(y = d$cases,        aff = d$AFF,
                e = d$expected)

lipInits <- list(b0 = 0, b1 = 0, tau = 1, sigma = 1)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                        data = lipData, inits = lipInits)


###### Compiling the model and MCMC ######
mcmc.out <- nimbleMCMC(code = lipCode, constants = lipConsts,
                       data = lipData, inits = lipInits,
                       nchains = 2, niter = 10000, summary = TRUE,
                       WAIC = TRUE, monitors = c("b0","b1","tau",
                                                 "sigma","p_b0","p_b1",
                                                 "p_tau","p_sigma","z0",
                                                 "z1","psi","p_psi"))

mcmc.out$summary

priorXposterior <- function(prior, posterior, chainn, niter=10000){
  df_beta <-  cbind(chainn[,posterior], rep(posterior,niter))
  df_pbeta <- cbind(chainn[, prior]   , rep(prior,niter))
  df_beta <-  rbind(df_beta,df_pbeta)
  
  df <- data.frame(chain = as.numeric(df_beta[,1]), beta = df_beta[,2])
  
  p <- ggplot(data = df, aes(x = chain, fill = beta, colour = beta))+
    geom_density(alpha = 0.4)+
    scale_x_continuous(expression(tau), 
                       expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_bw(base_size = 20)
  
  return(p)
}


priorXposterior("p_tau","tau",mcmc.out$samples$chain1)
priorXposterior("p_sigma","sigma",mcmc.out$samples$chain1)
priorXposterior("p_b0","b0",mcmc.out$samples$chain1)
priorXposterior("p_b1","b1",mcmc.out$samples$chain1)

table(mcmc.out$samples$chain1[,"z0"])
table(mcmc.out$samples$chain1[,"z1"])



tau <- 1
g <- 0.01
rcar_proper(n=1,mu = rep(0,53), adj = non.zero.adj$adj, num = non.zero.adj$num,
            tau = tau, gamma = g)
