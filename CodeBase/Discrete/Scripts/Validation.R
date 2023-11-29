library(data.table)
library(tscount)
library(sf)
library(tmap)
library(tidyverse)
library(gridExtra)

# Load data
setwd("~/Documents/MyResearch/CodeBase/Discrete/Scripts")
load("~/Documents/MyResearch/CodeBase/Discrete/Workspaces/data.RData")

# Prep
Z <- nrow(grid_map)
M <- 6
data_grid$id <- factor(data_grid$id, levels = 1:Z)
grid_map$id <- factor(grid_map$id, levels = 1:Z)

# Train and test split
T_END <- 2190
T_HORIZON <- 730
train_data <- data_grid[data_grid$time <= T_END, ]
test_data <- data_grid[data_grid$time > T_END & data_grid$time <= T_END + T_HORIZON, ]
mark_prop <- unclass(table(train_data$type))/nrow(train_data)

# TSGLM TRAIN
count_array = unclass(table(factor(train_data$time, levels = 1:T_END), train_data$id))
tsfit_pois_list <- rep(list(0), Z)
for(i in 1:Z){
  tsfit_pois_list[[i]] <- tsglm(ts(count_array[,i]), model = list(past_obs = 1, past_mean = 1),
                                  distr = "poisson")
}

# TSGLM PREDICT
baseline_array <- array(0, dim = c(T_HORIZON, Z))
count_array = unclass(table(factor(test_data$time, levels = (T_END+1):(T_END + T_HORIZON)), test_data$id))
for(j in 1:Z){
  newobs = ts(count_array[,j])
  baseline_array[,j] = predict(tsfit_pois_list[[j]], n.ahead = T_HORIZON, newobs = newobs, level = 0, global = TRUE)$pred
}

# MODEL PREDICT
load("~/Documents/MyResearch/CodeBase/Discrete/Workspaces/estimates.RData")
source("helpers.R")
#Parameters
mu_x_bw <- 4; mu_y_bw  <- 4
params <- list(mu_txym, alpha, g_t, g_xy, mu_x_breaks, mu_y_breaks, g_t_breaks, g_xy_breaks, t(gamma_cast),
               mu_x_bw, mu_y_bw, diff(g_t_breaks), diff(g_xy_breaks)) 
model_data <- data[, 1:4]
colnames(model_data) <- c('t', 'x', 'y', 'm')
truth_array  <- array(0, dim = c(T_HORIZON, Z, M))
model_array <- array(0, dim = c(T_HORIZON, Z, M))

for(i in (T_END+1):(T_END+T_HORIZON)){
  t_table <- unclass(table(data_grid$id[data_grid$time == i], data_grid$type[data_grid$time == i]))
  truth_array[i-T_END,,] = ifelse(is.null(ncol(table)), 0, t_table)
  history_i = model_data[model_data$t < i]
  sim_events = data.frame()
  for(j in 1:500){
    sim_events = rbind(sim_events, event_simulator(t = i, params = params, events = history_i))
  }
  sim_grid <- grid_map %>%
    st_join( st_as_sf(sim_events, coords = c('x', 'y'), crs = st_crs(grid_map)), left = TRUE) %>%
    arrange(t) %>%
    filter(is.finite(t)) 
  m_table <- unclass(table(sim_grid$id[sim_grid$t == i], sim_grid$m[sim_grid$t == i]))
  model_array[i-T_END,,] <- m_table/500
}

model_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], model_array[, z, m], log = T))))
baseline_array[baseline_array < 1e-16] <- 1e-16
baseline_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], baseline_array[, z]*mark_prop[m], log = T))))
poisson_rate_3d = unclass(table(train_data$id, train_data$type))/T_END
poisson_rate_3d[poisson_rate_3d < 1e-16] <- 1e-16
poisson_scores_3d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_3d[z, m], log = T))))
poisson_rate_2d = as.numeric(table(train_data$id))/(T_END*M)
poisson_scores_2d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d[z], log = T))))
poisson_rate_1d = nrow(train_data)/(T_END*M*Z)
poisson_scores_1d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_1d, log = T))))
poisson_rate_2d_clever = as.numeric(table(train_data$id))/T_END
poisson_scores_2d_clever = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d_clever[z]*mark_prop[m], log = T))))

sum(model_scores)
sum(baseline_scores)
sum(poisson_scores_2d_clever)
sum(poisson_scores_2d)
sum(poisson_scores_3d)
sum(poisson_scores_1d)

model_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                          function(z) scoring(response = truth_array[, z, m], pred = model_array[, z, m], 
                                                                                         distr = "poisson"))))

baseline_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                             function(z) scoring(response = truth_array[, z, m], pred = baseline_array[, z]*mark_prop[m], 
                                                                                            distr = "poisson"))))

poisson_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                            function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_2d_clever[z]*mark_prop[m], T_HORIZON), 
                                                                                           distr = "poisson"))))

pred_perf <- as.data.frame(cbind(rowSums(model_scores_array), rowSums(baseline_scores_array),rowSums(poisson_scores_array)))
colnames(pred_perf) <- c('DTMSTPP', 'TSGLM', 'POISSON')
pred_perf <- round(pred_perf, 4)
#pred_perf$scoring <- row.names(pred_perf)

pdf("~/Documents/MyResearch/CodeBase/Discrete/Figures/PredPerf.pdf", height=3, width=5)
grid.table(pred_perf[-c(5,6),])
dev.off()

# Attack probs plot
set.seed(1234)
ind_shp <- st_read("~/Documents/MyResearch/CodeBase/Discrete/Shapefiles/gadm40_IND_shp/gadm40_IND_0.shp",
                   stringsAsFactors = FALSE)
ind_shp <- st_simplify(ind_shp, preserveTopology = TRUE, dTolerance = 1000)

brekks <- c(0, 0.01, 0.05, 0.1, 1)
k <- 2899
  attacks <- st_as_sf(data[data$time == k, c('long', 'lat')], coords = c('long', 'lat'), crs = st_crs(ind_shp))
  grid_map$ATTACKS <- rowSums(1 - exp(-model_array[k-T_END, , ]))*10
  
  mapZH <- tm_shape(grid_map) +
    tm_polygons("ATTACKS", breaks = brekks, contrast=1, palette="Blues", 
                title = "Probability", labels = c('Less than 1%', '1% to 5%', '5% to 10%','10% and more'))
  
  pdf(paste0("~/Documents/MyResearch/CodeBase/Discrete/Figures/Validation/GridMapTruth", k, ".pdf"), height=6, width=8)
  print(qtm(grid_map, fill= "#b7c6e5") +
    tm_layout(scale = .5,
              legend.position = c(1.05,.775), 
              legend.title.size = 3,
              legend.height = 1,
              legend.text.size = 1.5,
              frame = F) + 
    mapZH + 
    tm_shape(attacks) + 
    tm_symbols(size = 1, alpha = 0.75, col = 'red', border.col = 'white'))
  dev.off()
