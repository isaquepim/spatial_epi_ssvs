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

ggplot()+
  geom_sf(data = ebola, aes(fill = include))+
  scale_fill_manual(values=c(OUT = "#74736E", IN = "#FAAC5A"))+
  labs(fill="Incluir")+
  theme_void()

ggplot()+
  geom_sf(data = ebola.raw, aes(fill = SIR))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "SIR")+
  theme_void()

ggplot()+
  geom_sf(data = ebola.raw, aes(fill = Case_count))+
  scale_fill_gradient(low = "#FEEDDE", high = "#A63603", name = "Casos")+
  theme_void()

#####################
#neighbour list



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


#nimble was not accepting ebola_IN.covariates, so had to find a quick fix
cov.ebola <- matrix(rep(0,N.countries*N.covariates),
                    nrow = N.countries, ncol = N.covariates)

for(i in 1:N.countries){
  for(j in 1:N.covariates){
    cov.ebola[i,j] <- as.numeric(ebola_IN.covariates[i,j])
  }
}



