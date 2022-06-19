library(ggplot2)
library(patchwork)
library(sf)


plot_cart <- function(map, limits){
  void <- ggplot() + theme_void()
  
  prior1 <- ggplot() + 
    geom_sf(data = map, aes(fill = prior.lower)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  prior2 <- ggplot() + 
    geom_sf(data = map, aes(fill = prior.mean)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  prior3 <- ggplot() + 
    geom_sf(data = map, aes(fill = prior.upper)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  reference <- ggplot() + 
    geom_sf(data = map, aes(fill = SIR)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  post1 <- ggplot() + 
    geom_sf(data = map, aes(fill = posterior.lower)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  post2 <- ggplot() + 
    geom_sf(data = map, aes(fill = posterior.mean)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  post3 <- ggplot() + 
    geom_sf(data = map, aes(fill = posterior.upper)) + 
    theme_void() + 
    scale_fill_continuous(limits = limits,
                          oob = scales::squish)
  
  return((prior1|prior2|prior3)/
         (void|reference|void)/
         (post1|post2|post3))
}








d$SIR <- d$cases/d$expected
map <- scotland$spatial.polygon
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)

map <- st_as_sf(map)


### Set CRS to WGS84
map <- st_transform(map, crs = 4326)
void <- ggplot() + theme_void()
p <- 
  ggplot() + 
  geom_sf(data = map, aes(fill = SIR)) + 
  theme_void() + 
  scale_fill_continuous(limits = c(0,5),
                        oob = scales::squish)
(p | p | p) / (void|p|void) / (p | p | p)
