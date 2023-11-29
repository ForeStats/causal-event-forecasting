library(data.table)
library(tscount)
library(sf)
library(tmap)
library(tidyverse)
library(gridExtra)
library(PRROC)
library(hrbrthemes)
library(Rfast)
library(ggpubr)

# Load data
setwd("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts")
load("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/data.RData")

# Prep
Z            <- nrow(grid_map)
M            <- nrow(type_meta_data)
data_grid$id <- factor(data_grid$id, levels = 1:Z)
grid_map$id  <- factor(grid_map$id, levels = 1:Z)

# Train and test split
T_END      <- 8*52*7
T_HORIZON  <- max(data$time) - T_END
train_data <- data_grid[data_grid$time <= T_END, ]
test_data  <- data_grid[data_grid$time > T_END & data_grid$time <= T_END + T_HORIZON, ]
mark_prop  <- as.numeric(table(train_data$type))/nrow(train_data)
zone_prop  <- as.numeric(table(train_data$id))/nrow(train_data)
zone_prop[zone_prop < 1e-16] <- 1e-16
zone_prop  <- zone_prop/sum(zone_prop)

# TSGLM TRAIN 1D
count_array    <- as.numeric(table(factor(train_data$time, levels = 1:(T_END) )))
tsfit_pois_1d  <- tsglm(ts(count_array), model = list(past_obs = 1:2, past_mean = 1:4), distr = "poisson")
tsfit_pois_1d$logLik

# TSGLM TRAIN 2D
count_array   <- unclass(table(factor(train_data$time, levels = 1:(T_END)), factor(train_data$type, levels = 1:M)))
tsfit_pois_2d <- apply(count_array, 2, FUN = function(z) tsglm(ts(z), model = list(past_obs = 1:2, past_mean = 1:4),
                                               distr = "poisson"))

# TSGLM PREDICT 2D
count_array    <- unclass(table(factor(test_data$time, levels = (T_END + 1):(T_END + T_HORIZON) ), factor(test_data$type, levels = 1:M)))
baseline_array <- sapply(1:5, FUN = function(z) as.numeric(predict(tsfit_pois_2d[[z]], n.ahead = T_HORIZON, 
                                                                   newobs = ts(count_array[,z]), level = 0, global = TRUE)$pred) )

# save data
save.image(file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/L0/day/workspaceL0.Rdata')

# MODEL PREDICT
load("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/L0/estimatesL0.RData")
source("helpers.R")
#Parameters
params       <- list(mu_tm, alpha, g_t, g_t_breaks, t(gamma_cast), diff(g_t_breaks)) 
model_data   <- data[, c(3, 6)]
colnames(model_data) <- c('t', 'm')

truth_array  <- unclass(table( factor(test_data$time, levels = (T_END + 1):(T_END + T_HORIZON)), 
                               test_data$id, factor(test_data$type, levels = 1:M) ))

N_sim <- 500
sim_zk_list = lapply((T_END + 1):(T_END + T_HORIZON), FUN = function(z){
  history_i = model_data[model_data$t <= (z-1)]
  
  print(z)
  sim_z_list = lapply(1:N_sim, FUN = function(k){
    sim_k <- event_simulator_L0(start_t = (z-1), t_width = 1, params = params, events = history_i)
    if(nrow(sim_k) > 0){
      return(cbind(sim_id = k, sim_k))
    }
  })
  return(cbind(bin_id = z, Reduce(rbind, sim_z_list)))
})
sim_events <- Reduce(rbind, sim_zk_list)
sim_grid <- sim_events

model_array <- unclass(table( factor(sim_grid$bin_id, levels = (T_END + 1):(T_END + T_HORIZON)), 
                              factor(sim_grid$m, levels = 1:M) ))
model_array <- model_array/N_sim

poisson_array <- as.numeric(table(train_data$type))/(T_END)
poisson_array <- matrix(rep(poisson_array, T_HORIZON), nrow = T_HORIZON, ncol = M, byrow = TRUE)

save.image(file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/L0/day/workspaceL0.Rdata')

# Validation

# 2d models
model_scores    = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], model_array[,m]*zone_prop[z], log = T))))
baseline_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], baseline_array[,m]*zone_prop[z], log = T))))
poisson_scores  = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_array[,m]*zone_prop[z], log = T))))

sum(model_scores)
sum(baseline_scores)
sum(poisson_scores)

model_scores_array    = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                               function(z) colSums(scoring(response = truth_array[, z, m], pred = model_array[,m]*zone_prop[z], 
                                                     distr = "poisson", individual = TRUE)))) )

baseline_scores_array = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                               function(z) colSums(scoring(response = truth_array[, z, m], pred = baseline_array[,m]*zone_prop[z], 
                                                     distr = "poisson", individual = TRUE)))) )

poisson_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                               function(z) colSums(scoring(response = truth_array[, z, m], pred = poisson_array[,m]*zone_prop[z], 
                                                     distr = "poisson", individual = TRUE)))) )

pred_perf <- as.data.frame(cbind(rowSums(model_scores_array), rowSums(baseline_scores_array), rowSums(poisson_scores_array)))
colnames(pred_perf) <- c('DTMSTPP', 'TSGLM', 'POISSON')
pred_perf <- round(pred_perf/(T_HORIZON), 4)

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/L0/day/PredPerfL0.pdf", height=3, width=5)
grid.table(pred_perf[-c(5,6), -2])
dev.off()

# START SKIP
# 1d models
model_array <- as.numeric(table( factor(sim_grid$bin_id, levels = (T_END + 1):(T_END + T_HORIZON)) ))
model_array <- model_array/N_sim

# TSGLM PREDICT 1D
count_array    <- as.numeric(table(factor(test_data$week, levels = (T_END + 1):( (T_END + T_HORIZON) ) )))
baseline_array <- as.numeric(predict(tsfit_pois_1d, n.ahead = T_HORIZON, newobs = ts(count_array), level = 0, global = TRUE)$pred)

poisson_array  <- rep(nrow(train_data)/(T_END), T_HORIZON)

model_scores    = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], model_array*mark_prop[m]*zone_prop[z], log = T))))
baseline_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], baseline_array*mark_prop[m]*zone_prop[z], log = T))))
poisson_scores  = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_array*mark_prop[m]*zone_prop[z], log = T))))

sum(model_scores)
sum(baseline_scores)
sum(poisson_scores)

hist((model_array %*% t(mark_prop) * max(zone_prop)))
hist((baseline_array %*% t(mark_prop) * max(zone_prop)))
hist((poisson_array %*% t(mark_prop) * max(zone_prop)))
max(diff(model_array) %*% t(mark_prop) * max(zone_prop))
# END SKIP

# START PAI 
hist((model_array * max(zone_prop)), n = 30)
hist((baseline_array * max(zone_prop)), n = 30)
hist((poisson_array * max(zone_prop)), n = 30)
hist(diff(model_array) * max(zone_prop), n = 100)

# Prediction Accuracy Index
N_PAI <- 100
model_pais <- baseline_pais <- poisson_pais <- array(0, dim = c(N_PAI, N_PAI, 2))
theta <- (seq(0, 1, 1/N_PAI)^3)*0.04
zeta <- (seq(0, 1, 1/N_PAI)^3)*0.02
for(i in 1:N_PAI){
  for(j in 1:N_PAI){
    model_pais[i,j,]    <- calculate_PAI_2d(truth_array = truth_array, prob_array = model_array,
                                            zone_prop = zone_prop, theta = theta[i], zeta = zeta[j])
    baseline_pais[i,j,] <- calculate_PAI_2d(truth_array = truth_array, prob_array = baseline_array,
                                            zone_prop = zone_prop, theta = theta[i], zeta = zeta[j])
    poisson_pais[i,j,]  <- calculate_PAI_2d(truth_array = truth_array, prob_array = poisson_array,
                                            zone_prop = zone_prop, theta = theta[i], zeta = zeta[j])
  }
  print(i)
}

pai_list  <- list(DTMSTPP = model_pais, TSGLM = baseline_pais, POISSON = poisson_pais)
plot_data <- do.call(rbind, lapply( 1:3, FUN = function(z){
  
  pai_mat <- matrix(pai_list[[z]], N_PAI*N_PAI, 2)
  pai_mat[,1] <- pai_mat[,1]/sum(truth_array)
  pai_mat[,2] <- pai_mat[,2]/(T_HORIZON*Z*M)
  pai_mat <- cbind(pai_mat, pai_mat[,1]/pai_mat[,2])
  pais_mat <- data.table(pai_mat)
  colnames(pais_mat) <- c('hit_rate', 'area_perc', 'PAI')
  pais_mat <- pais_mat[, .SD[which.max(PAI)], by = c('hit_rate','area_perc')]
  pais_mat$hit_rate_bin <- cut(pais_mat$hit_rate, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
  #pais_mat$area_perc_bin <- cut(pais_mat$area_perc, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
  pais_mat <- pais_mat[, .SD[which.max(PAI)], by = c('hit_rate_bin')]
  
  return(cbind(pais_mat[, c('hit_rate', 'area_perc', 'PAI')], names(pai_list)[z] ))
}))

# plot data
colnames(plot_data) <- c('HIT.RATE', 'AREA.PERC', 'PAI', 'MODEL')
setorder(plot_data, MODEL, AREA.PERC)

p1 <- ggplot(plot_data[plot_data$MODEL != 'TSGLM',], aes(x=AREA.PERC, y=HIT.RATE, colour = MODEL)) +
  geom_line() + 
  xlab("% AREA") +
  ylab("HIT RATE") +
  xlim(0, 0.5) +
  ylim(0.2, 1) +
  theme_ipsum(base_size = 16,
              axis_title_size = 16) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom')

p2 <- ggplot(plot_data[plot_data$MODEL != 'TSGLM',], aes(x=HIT.RATE, y=PAI, colour = MODEL)) +
  geom_line() + 
  xlab("HIT RATE") +
  ylab("PAI") +
  xlim(0.2, 1) +
  ylim(0, 20) +
  theme_ipsum(base_size = 16,
              axis_title_size = 16) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom')

ggarrange(p1, p2)

save.image(file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/L0/day/workspaceL0.Rdata')

# END PAI

set.seed(123)

# Attack probs plot
#for(k in 1:T_HORIZON){
for(k in 9){
  attacks <- st_as_sf(data[data$time == (T_END + k) & data$type == 3, c('longitude', 'latitude')], 
                      coords = c('longitude', 'latitude'), crs = st_crs(ind_shp))
  
  if(nrow(attacks) > 0){

    prob_vec <- 1 - exp(-model_array[k, 3]*zone_prop)
    grid_map$ATTACKS <- prob_vec/max(prob_vec)
    
    mapZH <- tm_shape(grid_map) +
      tm_polygons("ATTACKS", style = 'pretty', contrast=1, palette="Blues", 
                  title = "Pr( Bombing )")
    
    pdf(paste0("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/L0/day/GridMapTruth", k, "_L0.pdf"), height=6, width=9)
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