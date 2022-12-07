library(nimble)
library(spdep)
library(SpatialEpi)
library(sp)
library(ggplot2)
library(INLA)
library(ggthemes)
library(RColorBrewer)
library(coda)
library(sf)

set.seed(123)
###### loading data ######

data(scotland)

d <- scotland$data
d$SIR <- d$cases/d$expected*100

###### setting up the model ######

UnstructCode <- nimbleCode({
  
  for (i in 1:N){
    # LogPoisson model
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + b1*aff[i] + u[i]
    u[i] ~ dnorm(0,tau.u)
    
    #area-specific SIR
    SIR[i] <- exp(b0 + b1*aff[i] + u[i])
    
    # residuals SIR and residual from observation
    resSIR[i] <- exp(u[i])
    res[i] <- (y[i] - mu[i])/sqrt(mu[i])
    
  }
  
  b0 ~ dnorm(0,0.001) #intercept prior
  b1 ~ dnorm(0,0.001) #aff coef. prior
  
  tau.u ~ dgamma(1,0.1) #unstructured effect prior
  var.u <- 1/sqrt(tau.u) #standard deviation for unstructured effect
  
  allSIR <- exp(b0) #Commom ground SIR
})

##Parameters

N.cases <- length(d$cases)

#mcmc

n.chains <- 4
n.burnin <- 2000
n.iter   <- 10000

monitors <- c("b0", "b1", "allSIR", "resSIR", "SIR",
              "res", "u", "var.u", "mu")

##Setup for NIMBLE model

UnstructConsts <- list(N = N.cases)
UnstructData <- list(y = d$cases, aff = d$AFF, e = d$expected)

# initiate chain from 4 differente places
UnstructInits <- list(
  list(b0 =  0, b1 = 0, u = rnorm(N.cases,0,1), tau.u = 1),   #chain 1
  list(b0 =  1, b1 = 1, u = rnorm(N.cases,0,1), tau.u = 1),   #chain 2
  list(b0 =  1, b1 = 5, u = rnorm(N.cases,0,1), tau.u = 1),   #chain 3
  list(b0 = -1, b1 = 5, u = rnorm(N.cases,0,1), tau.u = 1)   #chain 4
  
  )

###### Compiling and running the model ######




UnstrCodesamples <- nimbleMCMC(code = UnstructCode,
                               data = UnstructData,
                               constants = UnstructConsts, 
                               inits = UnstructInits,
                               monitors = monitors,
                               niter = n.iter,
                               nburnin = n.burnin,
                               nchains = n.chains, 
                               samplesAsCodaMCMC = TRUE, 
                               summary = TRUE, 
                               WAIC = TRUE)


#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(UnstrCodesamples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.05) 

UnstrCodesamples$summary$all.chains[c('b0', 'b1', 'var.u'),]




##################################
### Plotting posterior mean of SIR

map <- scotland$spatial.polygon
proj4string(map) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"

map <- spTransform(map, CRS("+proj=longlat +datum=WGS84 +no_defs"))

#shapefile to polygon dataframe

d$SIRunstr <- UnstrCodesamples$summary$all.chains[paste0("SIR[", 1:N.cases, "]"), "Median"]
br <- c(0,1,2,3,4,5)
br <- c(br, Inf)
d$dSIRunstr <- cut(d$SIRunstr, breaks = br, dig.lab = 5,
              include.lowest = TRUE)

rownames(d) <- d$county
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)
map <- st_as_sf(map)


mytheme <- theme(panel.grid.major = element_line(color = '#cccccc', 
                                                 linetype = 'dashed',
                                                 size = .3),
                 panel.background = element_rect(fill = 'aliceblue')
)

p1 <- ggplot()+
  geom_sf(data = map, aes(fill = SIRunstr), color = NA)+
  coord_sf(crs = 27700)+
  #scale_fill_colorblind(name = "SIR (%)")+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603",limits = c(0,5), name = "SIR")+
  mytheme
p1
