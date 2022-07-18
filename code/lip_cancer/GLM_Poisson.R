library(nimble)
library(SpatialEpi)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)

library(sf)

###### loading data ######

data(scotland)

d <- scotland$data
d$SIR <- d$cases/d$expected*100

###### setting up the model ######

lipPriors <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(5,1)
  
  for (i in 1:N){
    
    log(theta[i]) <- b0 + b1*aff[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
  
})

lipCode <- nimbleCode({
  
  b0 ~ dnorm(0,1)
  b1 ~ dnorm(5,1)
  
  for (i in 1:N){
    
    log(theta[i]) <- b0 + b1*aff[i]
    y[i] ~ dpois(e[i]*theta[i])
    
  }
})

lipConsts <- list(N = length(d$cases))
lipData <- list(y = d$cases, aff = d$AFF, e = d$expected)
lipInits <- list(b0 = 0, b1 = 0)
lipPriorData <- list(aff = d$AFF, e = d$expected)


lipModel <- nimbleModel(code = lipCode, name = "lip", constants = lipConsts,
                    data = lipData, inits = lipInits)

lipPriorModel <- nimbleModel(code = lipPriors, name = "priors", constants = lipConsts,
                             data = lipPriorData, inits = lipInits)


###### testing ######

set.seed(123)

###### Compiling the model and MCMC ######

# nimbleMCMC compiles the model. If u want to customize your MCMC
# compile it yourself using:
#C.LipModel <- compileNimble(lipModel)

monitors <- c("b0","b1")


cmodel  <- compileNimble(lipModel)
mcmc    <- buildMCMC(lipModel, monitors = monitors)
cmcmc   <- compileNimble(mcmc, project = lipModel)
samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)


prior.cmodel  <- compileNimble(lipPriorModel)
prior.mcmc    <- buildMCMC(lipPriorModel, monitors = monitors)
prior.cmcmc   <- compileNimble(prior.mcmc, project = lipPriorModel)
prior.samples <- runMCMC(prior.cmcmc, niter = 1000, nburnin = 500)



pp.plot(prior.samples = list(prior.samples),
        posterior.samples = list(samples))


########################
######SPATIAL PLOT######
########################


d$SIR <- d$cases/d$expected
n <- length(d$cases)

y_posterior <- predictive_check(monitors = monitors,
                                nSamp = nrow(samples),
                                n = length(d$cases),
                                cmodel = cmodel,
                                samples = samples)

y_prior <- predictive_check(monitors = monitors,
                            nSamp = nrow(prior.samples),
                            n = length(d$cases),
                            cmodel = prior.cmodel,
                            samples = prior.samples)


q05 <- function(x){quantile(x,0.05)}
q95 <- function(x){quantile(x,0.95)}

d$prior.mean <- colMeans(y_prior)/d$expected
d$prior.lower <- apply(y_prior,2,q05)/d$expected
d$prior.upper <- apply(y_prior,2,q95)/d$expected

d$posterior.mean <- colMeans(y_posterior)/d$expected
d$posterior.lower <- apply(y_posterior,2,q05)/d$expected
d$posterior.upper <- apply(y_posterior,2,q95)/d$expected


map <- scotland$spatial.polygon
proj4string(map) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"

map <- spTransform(map, CRS("+proj=longlat +datum=WGS84 +no_defs"))

#shapefile to polygon dataframe

br <- c(0,1,2,3,4,5)*100
br <- c(br, Inf)
d$dSIR <- cut(d$SIR, breaks = br, dig.lab = 5,
    include.lowest = TRUE)
d$dAFF <- cut(d$AFF*100,c(0,1,2,3,4,5)*5, include.lowest = TRUE)

rownames(d) <- d$county
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)
map <- st_as_sf(map)

plot_cart(map = map,limits = c(0,7))


#EPSG:27700
#OSGB 1936 / British National Grid --
#United Kingdom Ordnance Survey
mytheme <- theme(panel.grid.major = element_line(color = '#cccccc' 
                                                  ,linetype = 'dashed'
                                                  ,size = .3
                 )
                 ,panel.background = element_rect(fill = 'aliceblue')
)

lipMap <- ggplot()+
  geom_sf(data = map, aes(fill = dSIR), color = NA)+
  coord_sf(crs = 27700)+
  #scale_fill_colorblind(name = "SIR (%)")+
  scale_fill_brewer(palette = "Oranges", name = "SIR (%)")+
  mytheme
lipMap

lipMap <- ggplot()+
  geom_sf(data = map, aes(fill = dAFF), color = NA)+
  coord_sf(crs = 27700)+
  #scale_fill_colorblind(name = "SIR (%)")+
  scale_fill_brewer(palette = "Oranges", name = "AFF (%)")+
  mytheme
lipMap
