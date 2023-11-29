# Load packages
library(readxl)
library(data.table)
library(gridExtra)
library(ggplot2)
library(hrbrthemes)
library(sf)
library(tmap)
library(tidyverse)
library(RColorBrewer)
setwd("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts")
options(digits=10)

# Import data
rawdata <- read_excel("~/Documents/MyResearch/Datasets/Terrorism/globalterrorismdb_0221dist.xlsx")

# Filter Indian records and relevant columns  
colids <- c("eventid", "iyear", "imonth", "iday", "longitude", "latitude", "attacktype1")
data   <- data.table(rawdata[rawdata$country_txt == 'India', colids])

# Data cleaning
data$longitude[data$eventid == '201811250017'] <- 80.721442
data      <- data[is.finite(data$latitude) & is.finite(data$longitude), ]
data$date <- as.Date(paste(data$iyear, data$imonth, data$iday, sep = "/"))
data      <- data[is.finite(data$date)]
  
# Filter data from 4/1/2008 onwards
hist(data$date, "years")
hist(data$date, "months")
hist(data$date[data$date > '2008/03/31'], "months")
hist(data$date[data$date > '2008/03/31'], "weeks")
table(data$iyear, data$imonth)
table(data$iyear, data$attacktype1)
spp_data      <- data[data$date > '2008/03/31' & data$date < '2019/06/01']
spp_data$time <- as.numeric(spp_data$date - as.Date('2008/03/31'))
setorder(spp_data, date, eventid)
spp_data <- spp_data[,c(8:9, 5:7)]

# Set Mark types
spp_data$type_label <- 'OTHERS'
spp_data$type_label[spp_data$attacktype1 == 1] <- 'ASSASSINATION'
spp_data$type_label[spp_data$attacktype1 == 2] <- 'ARMED ASSAULT'
spp_data$type_label[spp_data$attacktype1 == 3] <- 'BOMBING/EXPLOSION'
spp_data$type_label[spp_data$attacktype1 == 6] <- 'KIDNAPPING'
spp_data$type_label[spp_data$attacktype1 == 7] <- 'FACILITY/INFRASTRUCTURE ATTACK'
spp_data <- spp_data[spp_data$type_label != 'OTHERS']
spp_data$type  <- factor(spp_data$type_label, labels = 1:5)
type_meta_data <- unique(spp_data[, c('type', 'type_label')])
setorder(type_meta_data, type)
spp_data <- spp_data[,c(1:4, 7, 6)]

snapdata <- head(spp_data[,-c(2,5)], 10)
snapdata[,2:3] <- round(snapdata[,2:3], 3)
colnames(snapdata)[4] <- 'type'
pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/DataSnap.pdf", height=4, width=7)
grid.table(snapdata)
dev.off()

spp_data$week  <- 1 + ((spp_data$time-1) %/% 7)
#hist(table(spp_data$week), n = 30)

# Filter columns
cols <- c('date', 'week', 'time', 'longitude', 'latitude', 'type')
data <- spp_data[,.SD, .SDcols = cols]
setkey(data, time)
data$type <- as.numeric(data$type)
rm(list = c("rawdata", "colids", "spp_data", "cols", "snapdata"))

# Map setup
set.seed(1234)
# ind_shp <- st_read("~/Documents/MyResearch/CodeBase/Shapefiles/India/gadm40_IND_0.shp",
#                    stringsAsFactors = FALSE)
ind_shp <- st_read("~/Documents/MyResearch/Projects/Terrorism/Shapefile/India_Boundary.shp",
                   stringsAsFactors = FALSE)
ind_shp <- st_simplify(ind_shp, preserveTopology = TRUE, dTolerance = 1000)

events   <- st_as_sf(as.data.frame(data)[,c('longitude', 'latitude')], coords = c('longitude', 'latitude'), crs = st_crs(ind_shp))
latlims  <- c(17.5, 26)
longlims <- c(79.25, 87.75)

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/DataMap.pdf", height=8, width=8)
ind_shp %>%
  ggplot() +
  geom_sf(lwd = 0.1, fill = NA) + 
  geom_sf(data = events, size = 0.5, alpha = 0.5, colour = 'lightblue') +
  coord_sf(expand = FALSE) +
  theme_void()
dev.off()

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/FilterMap.pdf", height=8, width=8)
ind_shp %>%
  ggplot() +
  geom_sf(lwd = 0.1, fill = NA) + 
  geom_sf(data = events, size = 0.5, alpha = 0.5, colour = 'lightblue') +
  geom_rect(
    xmin = longlims[1],
    ymin = latlims[1],
    xmax = longlims[2],
    ymax = latlims[2],
    fill = NA, 
    colour = "black",
    size = 0.5
  ) +
  coord_sf(expand = FALSE) +
  theme_void()
dev.off()

# Filter Bihar Chattishgarh Square
data <- data[(data$latitude > latlims[1]) &
               (data$latitude < latlims[2]) &
               (data$longitude > longlims[1]) &
               (data$longitude < longlims[2]),]

my_bbox <- st_as_sfc(st_bbox(c(xmin = longlims[1],
                               ymin = latlims[1],
                               xmax = longlims[2],
                               ymax = latlims[2]), crs = st_crs(ind_shp)))

grid_map <- my_bbox %>%
  st_make_grid(cellsize = 0.5, offset = c(longlims[1], latlims[1])) %>%
  st_sf() %>%
  mutate(id = row_number()) %>%
  cbind(do.call(rbind, st_make_grid(my_bbox, cellsize = 0.5, what = "centers", offset = c(longlims[1], latlims[1])))) %>%
  mutate(id = row_number())

tmap_format_add(frame = F, name = "frameless")
qtm(grid_map$geometry, fill = "white", format = "frameless", borders =  'black')

## Plotting
by_month <- function(x,n=1){
  seq(min(x,na.rm=T),max(x,na.rm=T),by=paste0(n," months"))
}

ggplot(data, aes(x = date)) +
  geom_histogram(breaks = by_month(data$date)) + 
  xlab("month") + 
  ylab("") + 
  ggtitle("Histogram of events per month") +
  theme_ipsum_rc(base_family = "Roboto Condensed",
                 base_size = 16,
                 axis_title_size = 16)

# Seasonality columns
data$year  <- year(data$date)
data$month <- month(data$date)

count_data      <- dcast.data.table(data, formula = month ~ year, fun = length)
count_data_plot <- melt(count_data, 1)
count_data_plot <- count_data_plot[count_data_plot$value > 0]

ggplot(count_data_plot,  aes(x=factor(month), y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  xlab("month") + 
  ylab("") + 
  ylim(c(5, 55)) + 
  ggtitle("Boxplot of events per month") +
  theme(axis.text.x = element_blank()) + 
  theme_ipsum_rc(base_family = "Roboto Condensed",
                 base_size = 16,
                 axis_title_size = 16)

data$bint  <- 1 + ((data$time-1) %/% (7*52*6) )
#round(as.matrix(table(data$year, data$type))/c(table(data$year)),2)
round(as.matrix(table(data$bint, data$type))/c(table(data$bint)),2)

data$long <- data$longitude
data$lat  <- data$latitude
data$year <- NULL
data$month <- NULL

data_grid <- grid_map %>%
  st_join( st_as_sf(data, coords = c('longitude', 'latitude'), crs = st_crs(grid_map)), left = TRUE) %>%
  arrange(time) %>%
  filter(is.finite(time)) 

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/GridMap.pdf", height=6, width=6)
qtm(grid_map, fill= "#ffffff") +
  tm_layout(scale = .5,
            frame = F) + 
  tm_shape(st_as_sf(as.data.frame(data_grid)[,c('long', 'lat')], coords = c('long', 'lat'), crs = st_crs(ind_shp))) + 
  tm_symbols(size = 1, alpha = 0.75, col = 'lightblue', border.col = 'white')
dev.off()

data_grid_counts <- data_grid  %>%
  group_by(id) %>% 
  summarise(EVENTS = n())

quantile(data_grid_counts$EVENTS, seq(0, 1, 0.01))
brekks <- c(0, 1, 5, 20, 50, 200)

mapZH <- tm_shape(data_grid_counts) +
  tm_polygons("EVENTS", breaks = brekks, contrast=c(0.1, 1), palette="Blues", 
              title = "EVENTS", labels = c('1', '2 to 5', '6 to 20', '21 to 50', '51 and more'),
              interval.closure = "left")

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/GridMapCounts.pdf", height=6, width=9)
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
dt_events <- data$time
N      <- max(data$time)
p_hat  <- length(dt_events)/N
xlim   <- 100
observed_K <- vapply(1:xlim, function(z) InHomoDiscK(s = z, events = dt_events, N = N - xlim, p_hat = p_hat) - z, FUN.VALUE = 1.0)

# Simulations on non-clustered Poisson data
set.seed(1234)
Nsim <- 100

poisson_df <- data.frame(vapply(1:Nsim, function(k){
  ps_events <- sort(sample(1:N, size = length(dt_events), prob = rep(p_hat, N), replace = TRUE))
  vapply(1:xlim, function(z) InHomoDiscK(s = z, events = ps_events, N = N - xlim, p_hat = p_hat) - z, FUN.VALUE = 1.0)
}, FUN.VALUE = rep(1.0, xlim)))

myColours <- brewer.pal(3,"Dark2")

plot_data <- data.frame(x = seq(1, xlim, 1), Poisson= poisson_df)
plot_data <- melt(data.table(plot_data), id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Poisson'

#  Boxplot
p1 <- ggplot() + 
  geom_point(aes(x = factor(1:xlim), y = observed_K, fill = 'Observed'))  +
  scale_fill_manual(name = NULL, labels = c('Observed times'), values = 'black') +
  geom_boxplot(data = plot_data, aes(x = factor(x), y = value, colour = 'Poisson'), outlier.shape = NA) +
  scale_colour_manual(name = NULL,
                      labels = c('Poisson'),
                      values = c("Poisson"=myColours[1])) +
  scale_x_discrete(breaks = c(0, 25, 50 , 75)) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  t") +
  theme_bw() + 
  theme(legend.position="top",
        text=element_text(size=18, vjust = 0.9))

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/Kfuncs.pdf", width=10, height=8)
p1
dev.off()

data$long <- NULL
data$lat  <- NULL

save(data, type_meta_data, data_grid, grid_map, ind_shp, my_bbox, longlims, latlims, file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/data.Rdata')
