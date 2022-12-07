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

StructCode <- nimbleCode({
  
  for (i in 1:N.non_islands){
    # LogPoisson model with CAR errors
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(e[i]) + b0 + b1*aff[i] + s[i]
    
    #area-specific SIR
    SIR[i] <- exp(b0 + b1*aff[i] + s[i])
    
    # residuals SIR and residual from observation
    resSIR[i] <- exp(s[i])
    res[i] <- (y[i] - mu[i])/sqrt(mu[i])
    
  }
  
  for(i in 1:N.islands){
    # LogPoisson model with CAR errors
    y.i[i] ~ dpois(mu.i[i])
    log(mu.i[i]) <- log(e.i[i]) + b0 + b1*aff.i[i] + u.i[i]
    u.i[i] ~ dnorm(0,tau.u)
    
    #area-specific SIR
    SIR.i[i] <- exp(b0 + b1*aff.i[i] + u.i[i])
    
    # residuals SIR and residual from observation
    resSIR.i[i] <- exp(u.i[i])
    res.i[i] <- (y.i[i] - mu.i[i])/sqrt(mu.i[i])    
  }
  
  b0 ~ dnorm(0,0.001) #intercept prior
  b1 ~ dnorm(0,0.001) #aff coef. prior
  
  s[1:N.non_islands] ~ dcar_normal(adj[1:L], w[1:L], num[1:N.non_islands],
                                   tau.s, zero_mean = 1) # CAR
  

  
  tau.s ~ dgamma(1, 0.1) #CAR precision prior
  var.s <- 1/sqrt(tau.s) #standard deviation for structured effect
  
  tau.u ~ dgamma(1, 0.1) #islands precision prior
  var.u <- 1/sqrt(tau.u) #standard deviation for islands unstructured effect
  
  allSIR <- exp(b0) #Commom ground SIR
})

##Parameters

N.cases <- length(d$cases)

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
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
islands <- which(card(nb) == 0)

# separating islands form mainland
W.no_islands <- W[-islands, -islands]
cases.no_islands <- d$cases[-islands]
aff.no_islands <- d$AFF[-islands]
E.no_islands <- d$expected[-islands]

cases.islands <- d$cases[islands]
aff.islands <- d$AFF[islands]
E.islands <- d$expected[islands]

#matrix to adjacency list
LipAdj <- as.carAdjacency(W.no_islands)

N.no_islands <- length(cases.no_islands)
N.islands <- length(islands)

#mcmc

n.chains <- 4
n.burnin <- 2000
n.iter   <- 10000

monitors <- c("b0", "b1", "allSIR", "resSIR", "SIR", "SIR.i",
              "res", "s", "var.s", "mu")

##Setup for NIMBLE model

StructConsts <- list(N.islands = N.islands,
                       N.non_islands = N.no_islands,
                       L = length(LipAdj$adj),
                       adj = LipAdj$adj, w = LipAdj$weights, num = LipAdj$num)

StructData <- list(y.i = cases.islands, aff.i = aff.islands,
                     e.i = E.islands,
                     y = cases.no_islands, aff = aff.no_islands,
                     e = E.no_islands)

# initiate chain from 4 differente places
StructInits <- list(
  list(b0 =  0, b1 = 0, u.i = rnorm(N.islands,0,1),
       s = rnorm(N.no_islands,0,1), tau.u = 1,tau.s = 1),   #chain 1
  
  list(b0 =  1, b1 = 1, u.i = rnorm(N.islands,0,1),
       s = rnorm(N.no_islands,0,1), tau.u = 1,tau.s = 1),   #chain 2
  
  list(b0 =  1, b1 = 5, u.i = rnorm(N.islands,0,1),
       s = rnorm(N.no_islands,0,1), tau.u = 1,tau.s = 1),   #chain 3
  
  list(b0 = -1, b1 = 5, u.i = rnorm(N.islands,0,1),
       s = rnorm(N.no_islands,0,1), tau.u = 1,tau.s = 1)   #chain 4
  
)

###### Compiling and running the model ######




StrCodesamples <- nimbleMCMC(code = StructCode,
                               data = StructData,
                               constants = StructConsts, 
                               inits = StructInits,
                               monitors = monitors,
                               niter = n.iter,
                               nburnin = n.burnin,
                               nchains = n.chains, 
                               samplesAsCodaMCMC = TRUE, 
                               summary = TRUE, 
                               WAIC = TRUE)


#Rhat to attest convergence (Gelman-Rubin diagnostic)
GR.diag <- gelman.diag(StrCodesamples$samples, multivariate = FALSE)
#If TRUE, all parameters have Rhat less than 1.1,
#an indicative of convergence
all(GR.diag$psrf[,"Point est."] < 1.05) 

StrCodesamples$summary$all.chains[c('b0', 'b1', 'var.s'),]


##################################
### Plotting posterior mean of SIR

map <- scotland$spatial.polygon
d <- scotland$data

map.islands <- map[islands]
map <- map[-islands]
proj4string(map)<- proj4string(map.islands) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"

map <- spTransform(map, CRS("+proj=longlat +datum=WGS84 +no_defs"))
map.islands <- spTransform(map.islands, CRS("+proj=longlat +datum=WGS84 +no_defs"))

#shapefile to polygon dataframe
d.islands <- d[islands,]
d <- d[-islands,]

d$SIRposterior <- StrCodesamples$summary$all.chains[paste0("SIR[", 1:N.no_islands, "]"), "Median"]
d.islands$SIRposterior <- StrCodesamples$summary$all.chains[paste0("SIR.i[", 1:N.islands, "]"), "Median"]

rownames(d) <- d$county
rownames(d.islands) <- d.islands$county

map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)
map <- st_as_sf(map)

map.islands <- SpatialPolygonsDataFrame(map.islands, d.islands, match.ID = TRUE)
map.islands <- st_as_sf(map.islands)


mytheme <- theme(panel.grid.major = element_line(color = '#cccccc', 
                                                 linetype = 'dashed',
                                                 size = .3),
                 panel.background = element_rect(fill = 'aliceblue')
)

p2 <- ggplot()+
  geom_sf(data = map, aes(fill = SIRposterior), color = NA)+
  geom_sf(data = map.islands, aes(fill = SIRposterior), color = NA)+
  coord_sf(crs = 27700)+
  #scale_fill_colorblind(name = "SIR (%)")+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", limits = c(0,5), name = "SIR")+
  mytheme
p2
