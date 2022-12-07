library(ggplot2)
library(ggmap)
library(rgdal)
library(maps)
library(mapdata)
library(broom)
library(patchwork)
library(sf)


base_dir <- "D:\\Profile\\OneDrive - Imagem Geosistemas e Comercio LTDA\\Documentos\\FGV\\TCC\\spatial_epi_tcc\\code\\maps"


ES.shp <- readOGR(
  dsn = paste0(base_dir, '\\shapefiles'),
  layer = 'ES'
)

ES.fortified <- tidy(ES.shp)
ES.sf <- st_as_sf(ES.shp)

ggplot()+
  geom_sf(data = ES.sf)+
  theme_void()




n <- nrow(ES.sf)

g <- as.data.frame(matrix(nrow=(n*(n-1)/2),ncol=4))
colnames(g) <- c('from_x')
k <- 0
for(i in 1:(n-1)){
  for(j in (i+1):n){
    k <- k+1
    
  }
}
