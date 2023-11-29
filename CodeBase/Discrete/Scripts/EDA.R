# Load packages
library(readxl)
library(readr)
library(data.table)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(tmap)
library(tidyverse)
setwd("~/Documents/MyResearch/CodeBase/Discrete/Scripts")
options(digits=10)

# Import data
rawdata <- read_excel("~/Documents/MyResearch/CodeBase/Discrete/Data/globalterrorismdb_0221dist.xlsx")

# Filter Indian records and relevant columns
colids <- c("eventid", "iyear", "imonth", "iday", "latitude", "longitude", "attacktype1")
data <- data.table(rawdata[rawdata$country_txt == 'India', colids])

# Data cleaning
data$longitude[data$eventid == '201811250017'] <- 80.72144
data <- data[is.finite(data$latitude) & is.finite(data$longitude), ]
data[data$iday == 0]
data$date <- as.Date(paste(data$iyear, data$imonth, data$iday, sep = "/"))
data <- data[is.finite(data$date)]
  
# Filter data from 1/1/2012 onwards
#hist(data$date, "years")
#hist(data$date, "months")
#hist(data$date[data$date > '2011/12/31'], "months")
#hist(data$date[data$date > '2011/12/31'], "weeks")
table(data$iyear, data$imonth)
spp_data <- data[data$date > '2011/12/31']
spp_data$time <- as.numeric(spp_data$date - as.Date("2011/12/31"))
setorder(spp_data, date, eventid)
spp_data <- spp_data[,c(8:9, 5:7)]
rm(list = c("data", "rawdata", "colids"))

# De-duplication and type labels
spp_data <- spp_data[!duplicated(spp_data),]

spp_data$type_label <- 'OTHERS'
spp_data$type_label[spp_data$attacktype1 == 1] <- 'ASSASSINATION'
spp_data$type_label[spp_data$attacktype1 == 2] <- 'ARMED ASSAULT'
spp_data$type_label[spp_data$attacktype1 == 3] <- 'BOMBING/EXPLOSION'
spp_data$type_label[spp_data$attacktype1 == 6] <- 'KIDNAPPING'
spp_data$type_label[spp_data$attacktype1 == 7] <- 'FACILITY/INFRASTRUCTURE ATTACK'
spp_data$type <- factor(spp_data$type_label, labels = 1:6)
type_meta_data <- unique(spp_data[, c('type', 'type_label')])
spp_data <- spp_data[,c(1:4, 7, 6)]

latlims <- c(22, 26)
longlims <- c(83.5, 87.5)
data <- spp_data[(spp_data$latitude > latlims[1]) &
                   (spp_data$latitude < latlims[2]) &
                   (spp_data$longitude > longlims[1]) &
                   (spp_data$longitude < longlims[2]),]

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/DataSnap.pdf", height=4, width=7)
grid.table(head(data[,-c(2,5)], 10))
dev.off()

# Set unique lat/long
cols <- c('time', 'longitude', 'latitude', 'type')
data <- data[,.SD, .SDcols = cols]
setkey(data, time)
set.seed(12345)
data$longitude <- data$longitude + runif(nrow(data), min = -1e-6, max = 1e-6)
data$latitude <- data$latitude + runif(nrow(data), min = -1e-6, max = 1e-6)
length(unique(data$longitude))
length(unique(data$latitude))

data$type <- as.numeric(data$type)
setkey(type_meta_data, type)

# Map setup
set.seed(1234)
ind_shp <- st_read("~/Documents/MyResearch/CodeBase/Discrete/Shapefiles/gadm40_IND_shp/gadm40_IND_0.shp",
                   stringsAsFactors = FALSE)
ind_shp <- st_simplify(ind_shp, preserveTopology = TRUE, dTolerance = 1000)

events <- st_as_sf(as.data.frame(spp_data)[,c('longitude', 'latitude')], coords = c('longitude', 'latitude'), crs = st_crs(ind_shp))

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/FilterMap.pdf", height=8, width=8)
ind_shp %>%
  ggplot() +
  geom_sf(lwd = 0.1, fill = NA) + 
  geom_sf(data = events, size = 0.5, alpha = 0.5, colour = 'lightblue') +
  geom_rect(
    xmin = 83.5,
    ymin = 22,
    xmax = 87.5,
    ymax = 26,
    fill = NA, 
    colour = "black",
    size = 0.5
  ) +
  coord_sf(expand = FALSE) +
  theme_void()
dev.off()

my_bbox <- st_as_sfc(st_bbox(c(xmin = 83.5,
                     ymin = 22,
                     xmax = 87.5,
                     ymax = 26), crs = st_crs(ind_shp)))

grid_map <- my_bbox %>%
  st_make_grid(cellsize = 0.5, offset = c(83.5, 22)) %>%
  st_sf() %>%
  mutate(id = row_number()) %>%
  cbind(do.call(rbind, st_make_grid(my_bbox, cellsize = 0.5, what = "centers", offset = c(83.5, 22)))) %>%
  mutate(id = row_number())

tmap_format_add(frame = F, name = "frameless")
qtm(grid_map$geometry, fill = "white", format = "frameless", borders =  'black')

data$long <- data$longitude
data$lat <- data$latitude

data_grid <- grid_map %>%
  st_join( st_as_sf(data, coords = c('longitude', 'latitude'), crs = st_crs(grid_map)), left = TRUE) %>%
  arrange(time) %>%
  filter(is.finite(time)) 

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/GridMap.pdf", height=6, width=6)
qtm(grid_map, fill= "#ffffff") +
  tm_layout(scale = .5,
            frame = F) + 
  tm_shape(st_as_sf(as.data.frame(data_grid)[,c('long', 'lat')], coords = c('long', 'lat'), crs = st_crs(ind_shp))) + 
  tm_symbols(size = 1, alpha = 0.75, col = 'lightblue', border.col = 'white')
dev.off()

data_grid_counts <- data_grid  %>%
  group_by(id) %>% 
  summarise(ATTACKS = n())

quantile(data_grid_counts$ATTACKS, seq(0, 1, 0.01))
brekks <- c(0, 1, 4, 10, 30, 50)

mapZH <- tm_shape(data_grid_counts) +
  tm_polygons("ATTACKS", breaks = brekks, contrast=c(0.1, 1), palette="Blues", 
              title = "ATTACKS", labels = c('1', '2 to 4', '5 to 10', '11 to 30', '31 and more'),
              interval.closure = "left")

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/GridMapCounts.pdf", height=6, width=9)
qtm(grid_map, fill= "#ffffff") +
  tm_layout(scale = .5,
            legend.position = c(1.1,.75), 
            legend.title.size = 3,
            legend.height = 1.5,
            legend.text.size = 2,
            frame = F) + 
  mapZH
dev.off()

# Second order test for subregion
# probability of at-least one event per day
source("helpers.R")
events <- data$time
N <- max(data$time)
p_hat <- length(unique(events))/N
xlim <- 100
observed_K <- sapply(1:xlim, function(z) InHomoDiscK(s = z, events = unique(events), N = N - xlim, p_hat = p_hat) - z)

# Simulations on non-clustered Poisson data
set.seed(1234)
Nsim <- 100
poisson_df <- data.frame()

for(i in 1:Nsim){
  
  events <- sort(sample(1:N, size = length(unique(events)), prob = rep(p_hat, N)))
  poisson_values <- sapply(1:xlim, function(z) InHomoDiscK(s = z, events = events, N = N - xlim, p_hat = p_hat) - z)
  poisson_df <- rbind(poisson_df, poisson_values)
}

myColours <- brewer.pal(3,"Dark2")

plot_data = data.frame(x = seq(1, 100, 1), Poisson= t(poisson_df))
plot_data <- melt(data.table(plot_data), id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Poisson'

#  Boxplot
p1 <- ggplot() + 
  geom_point(aes(x = factor(1:100), y = observed_K, fill = 'Observed'))  +
  scale_fill_manual(name = NULL, labels = c('Observed times'), values = 'black') +
  geom_boxplot(data = plot_data, aes(x = factor(x), y = value, colour = 'Poisson'), outlier.shape = NA) +
  scale_colour_manual(name = NULL,
                      labels = c('Poisson'),
                      values = c("Poisson"=myColours[1])) +
  scale_x_discrete(breaks = c(0, 25, 50 , 75)) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  2t") +
  theme_bw() + 
  theme(legend.position="top",
        text=element_text(size=18, vjust = 0.9))

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/Kfuncs.pdf", width=10, height=8)
p1
dev.off()

save(data, type_meta_data, data_grid, grid_map, ind_shp, my_bbox, file = '~/Documents/MyResearch/CodeBase/Discrete/Workspaces/data.Rdata')
