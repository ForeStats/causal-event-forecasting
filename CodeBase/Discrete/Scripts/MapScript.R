# Load packages
library(sf)
library(tmap)
library(tidyverse)

# Load data
setwd("~/Documents/MyResearch/CodeBase/Discrete/Scripts")
load("~/Documents/MyResearch/CodeBase/Discrete/Workspaces/data.RData")

# Map setup
set.seed(1234)
ind_shp <- st_read("~/Documents/MyResearch/CodeBase/Discrete/Shapefiles/gadm40_IND_shp/gadm40_IND_0.shp",
                   stringsAsFactors = FALSE)
ind_shp <- st_simplify(ind_shp, preserveTopology = TRUE, dTolerance = 1000)
tmap_format_add(frame = F, name = "frameless")
qtm(ind_shp$geometry, fill = "#c2d8ba", format = "frameless", borders = 'white')
range(spp_data$latitude)
range(spp_data$longitude)

grid_map <- ind_shp %>%
  st_make_grid(cellsize = 0.2, offset = c(68.7, 8.1)) %>%
  st_sf() %>%
  mutate(id = row_number()) %>%
  cbind(do.call(rbind, st_make_grid(ind_shp, cellsize = 0.2, what = "centers", offset = c(68.7, 8.1)))) %>%
  st_intersection(ind_shp) %>%
  select(-c(4,5)) %>%
  mutate(id = row_number())

qtm(grid_map$geometry, fill = "#c2d8ba", format = "frameless", borders =  'white')

# grid_map$ncoord <- vapply(grid_map$geometry, FUN = function(z) 
#   ifelse(is.null(dim(unlist(z[[1]]))[1]), 0, dim(unlist(z[[1]]))[1]), FUN.VALUE = 1)
# grid_map$area <- as.numeric(st_area(grid_map$geometry))
# quantile(grid_map$area, seq(0, 0.3, 0.001))
# 
# grid_map$flag <- 0
# grid_map$flag[grid_map$ncoord == 5 & grid_map$area >= 4e8] <- 1
# 
# qtm(grid_map$geometry[grid_map$flag==1], fill = "#c2d8ba", format = "frameless", borders =  'white')

spp_data$long <- spp_data$longitude
spp_data$lat <- spp_data$latitude

data_grid <- grid_map %>%
  st_join( st_as_sf(spp_data, coords = c('longitude', 'latitude'), crs = st_crs(grid_map)), left = TRUE) %>%
  arrange(time) %>%
  filter(is.finite(time)) 

data_grid_counts <- data_grid  %>%
  group_by(id) %>% 
  summarise(ATTACKS = n())

quantile(data_grid_counts$ATTACKS, seq(0.8, 1, 0.01))
brekks <- c(0, 1, 4, 10, 50, 350)

mapZH <- tm_shape(data_grid_counts) +
  tm_polygons("ATTACKS", breaks = brekks, contrast=c(0.1, 1), palette="Blues", 
              title = "ATTACKS", labels = c('1', '2 to 4', '5 to 10', '11 to 50', '51 and more'),
              interval.closure = "left")

qtm(grid_map, fill= "#ffffff") +
  tm_layout(scale = .5,
            legend.position = c(0.65,.75), 
            legend.title.size = 3,
            legend.height = 1.5,
            legend.text.size = 2,
            frame = F) + 
  mapZH

qtm(ind_shp, fill= "#ffffff") +
  tm_layout(scale = .5,
            legend.position = c(0.85,.15), 
            legend.title.size = 3,
            legend.height = 1.5,
            legend.text.size = 2,
            frame = F) + 
  tm_shape(st_as_sf(as.data.frame(data_grid)[,c('long', 'lat')], coords = c('long', 'lat'), crs = st_crs(ind_shp))) + 
  tm_symbols(size = 0.5, alpha = 0.5, col = 'lightblue', border.col = 'white')
