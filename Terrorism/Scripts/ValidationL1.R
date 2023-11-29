library(data.table)
library(tscount)
library(sf)
library(tmap)
library(tidyverse)
library(gridExtra)
library(PRROC)
library(hrbrthemes)

# Load data
setwd("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts")
load("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/data.RData")

# Prep
Z            <- nrow(grid_map)
M            <- nrow(type_meta_data)
data_grid$id <- factor(data_grid$id, levels = 1:Z)
grid_map$id  <- factor(grid_map$id, levels = 1:Z)

# Train and test split
T_END      <- 56*52
T_HORIZON  <- max(data$time) - T_END
train_data <- data_grid[data_grid$time <= T_END, ]
test_data  <- data_grid[data_grid$time > T_END & data_grid$time <= T_END + T_HORIZON, ]
mark_prop  <- as.numeric(table(train_data$type))/nrow(train_data)
zone_prop  <- as.numeric(table(train_data$id))/nrow(train_data)
zone_prop[zone_prop < 1e-16] <- 1e-16
zone_prop  <- zone_prop/sum(zone_prop)

# TSGLM TRAIN 1D
count_array <- as.numeric(table(factor(train_data$time, levels = 1:T_END)))
tsfit_pois  <- tsglm(ts(count_array), model = list(past_obs = 1:4, past_mean = 1:4), distr = "poisson")
tsfit_pois$logLik

# TSGLM TRAIN 2D
# count_array = unclass(table(factor(train_data$time, levels = 1:T_END), train_data$id))
# tsfit_pois_list <- apply(count_array, 2, FUN = function(z) tsglm(ts(z), model = list(past_obs = 1, past_mean = 1),
#                                               distr = "poisson"))

# TSGLM PREDICT 1D
count_array    <- as.numeric(table(factor(test_data$time, levels = (T_END + 1):(T_END + T_HORIZON) )))
baseline_array <- as.numeric(predict(tsfit_pois, n.ahead = T_HORIZON, newobs = ts(count_array), level = 0, global = TRUE)$pred)

# save data
save.image(file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/workspace.Rdata')

# MODEL PREDICT
load("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/estimates.RData")
source("helpers.R")
#Parameters
params       <- list(mu_txym, alpha, g_t, g_xy, mu_x_breaks, mu_y_breaks, g_t_breaks, g_xy_breaks, t(gamma_cast),
                     diff(mu_x_breaks), diff(mu_y_breaks), diff(g_t_breaks), diff(g_xy_breaks)) 
model_data   <- data[, 3:6]
colnames(model_data) <- c('t', 'x', 'y', 'm')
model_array  <- array(0, dim = c(T_HORIZON, Z, M))

truth_array  <- unclass(table( factor(data_grid$time, levels = (T_END + 1):(T_END + T_HORIZON)), 
                          data_grid$id, factor(data_grid$type, levels = 1:M) ))

sim_zk_list = lapply((T_END + 1):(T_END + T_HORIZON), FUN = function(z){
  history_i = model_data[model_data$t < z]
  
  sim_z_list = lapply(1:1000, FUN = function(k){
    sim_k <- event_simulator(start_t = (z-1), t_width = 1, params = params, events = history_i)
    if(nrow(sim_k) > 0){
      return(cbind(sim_id = k, sim_k))
    }
  })
  return(cbind(time_id = z, Reduce(rbind, sim_z_list)))
})
sim_events <- Reduce(rbind, sim_zk_list)

sim_grid <- grid_map %>%
  st_join( st_as_sf(sim_events, coords = c('x', 'y'), crs = st_crs(grid_map)), left = TRUE) %>%
  arrange(t) %>%
  filter(is.finite(t)) 

hist(table(sim_grid$time_id)/1000, n = 30)
hist(baseline_array, n = 30)

rbind(table(train_data$type)/T_END, table(test_data$type)/T_HORIZON)

rbind(table(sim_grid$m)/1000, table(test_data$type), sum(baseline_array)*mark_prop)
rbind(table(sim_grid$id)/1000, table(test_data$id), sum(baseline_array)*zone_prop)

model_array <- unclass(table( factor(sim_grid$time_id, levels = (T_END + 1):(T_END + T_HORIZON)), 
                               sim_grid$id, factor(sim_grid$m, levels = 1:M) ))
model_array <- model_array/1000

length(which(model_array > 0))
model_array_org <- model_array
save.image(file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/workspace.Rdata')

model_array <- model_array_org
model_array[model_array < 1e-16] <- 1e-16
model_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], model_array[, z, m], log = T))))
length(which(baseline_array > 1e-16))
baseline_array[baseline_array < 1e-16] <- 1e-16
baseline_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], baseline_array*mark_prop[m]*zone_prop[z], log = T))))
model_rate_clever = as.numeric(table(factor(sim_grid$time_id, levels = (T_END + 1):(T_END + T_HORIZON))))/1000
model_clever_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], model_rate_clever*mark_prop[m]*zone_prop[z], log = T))))
poisson_rate_3d = unclass(table(train_data$id, train_data$type))/T_END
poisson_rate_3d[poisson_rate_3d < 1e-16] <- 1e-16
poisson_scores_3d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_3d[z, m], log = T))))
poisson_rate_2d = as.numeric(table(train_data$id))/(T_END*M)
poisson_rate_2d[poisson_rate_2d < 1e-16] <- 1e-16
poisson_scores_2d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d[z], log = T))))
poisson_rate_1d = nrow(train_data)/(T_END*M*Z)
poisson_scores_1d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_1d, log = T))))
poisson_rate_2d_clever = as.numeric(table(train_data$id))/T_END
poisson_rate_2d_clever[poisson_rate_2d_clever < 1e-16] <- 1e-16
poisson_scores_2d_clever = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d_clever[z]*mark_prop[m], log = T))))

sum(model_scores)
sum(model_clever_scores)
sum(baseline_scores)
sum(poisson_scores_2d_clever)
sum(poisson_scores_2d)
sum(poisson_scores_3d)
sum(poisson_scores_1d)

model_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                          function(z) scoring(response = truth_array[, z, m], pred = model_array[, z, m], 
                                                                                         distr = "poisson"))))

baseline_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                             function(z) scoring(response = truth_array[, z, m], pred = baseline_array*mark_prop[m]*zone_prop[z], 
                                                                                            distr = "poisson"))))

poisson_scores_array_1d  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                                  function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_1d, T_HORIZON), 
                                                                                              distr = "poisson"))))

poisson_scores_array_2d  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                                    function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_2d[z], T_HORIZON), 
                                                                                              distr = "poisson"))))

poisson_scores_array_3d  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                            function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_3d[z, m], T_HORIZON), 
                                                                                           distr = "poisson"))))

poisson_scores_array_2d_clever  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                             function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_2d_clever[z]*mark_prop[m], T_HORIZON), 
                                                                                           distr = "poisson"))))

pred_perf <- as.data.frame(cbind(rowSums(model_scores_array), rowSums(baseline_scores_array), rowSums(poisson_scores_array_1d),
                                 rowSums(poisson_scores_array_2d), rowSums(poisson_scores_array_3d), rowSums(poisson_scores_array_2d_clever)))
colnames(pred_perf) <- c('DTMSTPP', 'TSGLM', 'POISSON_1D', 'POISSON_2D', 'POISSON_3D', 'POISSON_2D_CLEVER')
pred_perf <- round(pred_perf, 4)

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/PredPerf.pdf", height=3, width=5)
grid.table(pred_perf[-c(5,6), c(1,2,5)])
dev.off()

# Validation Classification
test_data$week <- factor(test_data$week, levels = 313:416)
truth_cast <- dcast.data.table(data.table(test_data), formula = week + id ~ type, value.var = 'time', 
                             fun.aggregate = length, drop = FALSE)
colnames(truth_cast)[3:8] <- paste0('M_', 1:6) 
truth_cast <- truth_cast[, lapply(.SD, FUN = function(z) ifelse(z > 0, 1, 0)), by = c('week', 'id'), .SDcols = c(paste0('M_', 1:6))]

# TSGLM PREDICT
baseline_array <- array(0, dim = c(T_HORIZON, Z))
count_array = unclass(table(factor(test_data$week, levels = (T_END + 1):(T_END + T_HORIZON)), test_data$id))
for(j in 1:Z){
  newobs = ts(count_array[,j])
  baseline_array[,j] = predict(tsfit_pois_list[[j]], n.ahead = T_HORIZON, newobs = newobs, level = 0, global = TRUE)$pred
}
baseline_cast <- 1 - exp(-baseline_array)
colnames(baseline_cast) <- 1:Z
baseline_cast <- data.table(cbind(week = 1:104, baseline_cast))
baseline_melt <- melt.data.table( baseline_cast, 1, variable.name = 'id', value.name = 'prob')
setorder(baseline_melt, week, id)

# Poisson PREDICT
poisson_rate_3d = unclass(table(train_data$id, train_data$type))/T_END
colnames(poisson_rate_3d) <- paste0('M_', 1:6)
poisson_rate_3d <- 1 - exp(-poisson_rate_3d)
poisson_rate_3d <- cbind(id = 1:289, poisson_rate_3d)
poisson_melt <- data.table(cbind(week = rep(313:416, each = 289), 
                                 poisson_rate_3d[rep(seq_len(nrow(poisson_rate_3d)), times = 104), ]))
setorder(poisson_melt, week, id)

sim_grid <- grid_map %>%
  st_join( st_as_sf(sim_events, coords = c('x', 'y'), crs = st_crs(grid_map)), left = TRUE) %>%
  arrange(w, sim_id, t) %>%
  filter(is.finite(t))

sim_grid$w <- factor(sim_grid$w, levels = 313:416)
sim_grid$sim_id <- factor(sim_grid$sim_id, levels = 1:500)
sim_grid$m <- factor(sim_grid$m, levels = 1:M)

sim_cast <- dcast.data.table(data.table(sim_grid), formula = w + id + sim_id ~ m, value.var = 't', 
                             fun.aggregate = length, drop = FALSE)
colnames(sim_cast)[4:9] <- paste0('M_', 1:6) 
sim_cast_agg <- sim_cast[, lapply(.SD, FUN = function(z) length(which(z>0))/500), by = c('w', 'id'), .SDcols = c(paste0('M_', 1:6))]
mstpp_preds <- sim_cast_agg
colnames(mstpp_preds)[1] <- 'week'

# ROC and PR Curves
df.plot <- data.frame(
                      Truth = melt(truth_cast, id.vars = 1:2)[,4], 
                      MSTPP = melt(mstpp_preds, id.vars = 1:2)[,4],
                      POISSON = melt(poisson_melt, id.vars = 1:2)[,4])
colnames(df.plot) <- c('Truth', 'MSTPP', 'POISSON')

ggprcs(prcs = list(MSTPP = pr.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$MSTPP,  curve = TRUE),
                   #TSGLM = pr.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$TSGLM,  curve = TRUE),
                   POISSON = pr.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$POISSON,  curve = TRUE)))

ggrocs(rocs = list(MSTPP = roc.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$MSTPP,  curve = TRUE),
                   #TSGLM = roc.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$TSGLM,  curve = TRUE),
                   POISSON = roc.curve(weights.class0 =  df.plot$Truth, scores.class0 = df.plot$POISSON,  curve = TRUE)))

k <- 3
df.plot <- data.frame(truth_cast[, c(1,2)], 
                      Truth = data.frame(truth_cast)[, k+2], 
                      MSTPP = data.frame(mstpp_preds)[, k+2])
df.plot$MSTPP <- (df.plot$MSTPP - min(df.plot$MSTPP))/(max(df.plot$MSTPP) - min(df.plot$MSTPP))

# Attack probs plot
for(k in 313:416){
  attacks <- st_as_sf(data[data$week == k & data$type == 3, c('longitude', 'latitude')], 
                      coords = c('longitude', 'latitude'), crs = st_crs(ind_shp))
  
  if(nrow(attacks) > 0){
    grid_map$ATTACKS <- df.plot$MSTPP[df.plot$week == k]
    
    mapZH <- tm_shape(grid_map) +
      tm_polygons("ATTACKS", style = 'pretty', contrast=1, palette="Blues", 
                  title = "Pr( Bombing )")
    
    pdf(paste0("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/Validation/GridMapTruth", k, ".pdf"), height=6, width=9)
    print(qtm(grid_map, fill= "#b7c6e5") +
            tm_layout(scale = .5,
                      legend.position = c(1.05,.775), 
                      legend.title.size = 3,
                      legend.height = 1,
                      legend.text.size = 1.5,
                      frame = F) + 
            mapZH + 
            tm_shape(attacks) + 
            tm_symbols(size = 1, alpha = 0.9, col = 'red', border.col = 'white'))
    dev.off()
  }
}