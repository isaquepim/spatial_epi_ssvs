d <- scotland$data[,c("county.names", "cases", "expected", "AFF")]
names(d) <- c("county", "Y", "E", "AFF")
d$SIR <- d$Y / d$E
d$prior <- media / d$E

sapply(slot(map, "polygons"), function(x){slot(x, "ID")})

rownames(d) <- d$county
map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)


library(leaflet)
library(rgeos)
l <- leaflet(map) %>% addTiles()

pal <- colorNumeric(palette = "YlOrRd", domain = map$SIR)

l %>%
  addPolygons(
    color = "grey", weight = 1,
    fillColor = ~ pal(SIR), fillOpacity = 0.5
  ) %>%
  addLegend(
    pal = pal, values = ~SIR, opacity = 0.5,
    title = "SIR", position = "bottomright"
  )

l <- leaflet(map) %>% addTiles()

pal <- colorNumeric(palette = "YlOrRd", domain = map$prior)

l %>%
  addPolygons(
    color = "grey", weight = 1,
    fillColor = ~ pal(prior), fillOpacity = 0.5
  ) %>%
  addLegend(
    pal = pal, values = ~prior, opacity = 0.5,
    title = "prior", position = "bottomright"
  )