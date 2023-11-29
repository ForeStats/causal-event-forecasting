# Load packages
library(Hmisc)
library(distances)
library(data.table)
library(ggplot2)
library(ggpubr)

# Load data
setwd("~/Documents/MyResearch/CodeBase/Discrete/Scripts")
load("~/Documents/MyResearch/CodeBase/Discrete/Workspaces/data.RData")

# times and x,y
T_END <- 730
data <- data.table(spp_data[spp_data$time <= T_END, ])
N <- nrow(data)
data <- cbind(index = 1:N, data)
row.names(data) <- NULL

# Initialise the P matrix
p_mat <- matrix(1, N, N)
p_mat[upper.tri(p_mat)] <- 0
p_mat <- p_mat/(1:N)

# Algorithm
# Background bandwidth
range(data$longitude)
range(data$latitude)
mu_t_bw <- T_END
mu_x_bw <- 30; mu_y_bw  <- 30
mu_x_breaks <- seq(68, 98, mu_x_bw)
mu_y_breaks <- seq(8, 38, mu_y_bw)

data$bint <- 1
data$binx <- as.numeric(cut2(data$longitude, cuts = mu_x_breaks, digits = 7))
data$biny <- as.numeric(cut2(data$latitude, cuts = mu_y_breaks, digits = 7))
#mu_df <- data[, .N, by = c('bint', 'binx', 'biny')]
#setorder(mu_df, bint, binx, biny)

# Excitation bandwidth
g_t_dist_mat <- as.matrix(distances(data$time))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

g_t_breaks = unique(c(0, 1, 4, seq(7, 49, 7)))
g_t_bw = diff(g_t_breaks)

data_g_t_intervals = apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0

g_xy_dist_mat <- as.matrix(distances(data[,c('longitude', 'latitude')]))
g_xy_dist_mat[upper.tri(g_xy_dist_mat)] <- NA

g_xy_breaks = c(0, 0.001, 0.01, 0.1)
g_xy_bw = diff(g_xy_breaks)

data_g_xy_intervals = apply(g_xy_dist_mat, 2, FUN = cut, breaks = g_xy_breaks, labels = F, include.lowest = TRUE)
data_g_xy_intervals[data_g_t_intervals == 0] <- 0
data_g_xy_intervals[is.na(data_g_xy_intervals)] <- 0
table(data_g_xy_intervals)

delta = 1
counter = 1
mu_data <- data.table(data[,c('index', 'bint', 'binx', 'biny')], p_ii = diag(p_mat))

while(delta > 1e-3){
  # Update mu
  mu_vec <- diag(p_mat)
  mu_data$p_ii <- mu_vec
  mu_txy = mu_data[, list(mu = sum(p_ii)/(mu_t_bw*mu_x_bw*mu_y_bw)), by = c('bint', 'binx', 'biny')]
    
  # Update triggering function
  offspring_sum <- N - sum(mu_vec)
  alpha <- offspring_sum/N
  
  g_df <- data.table(t_intvl = as.vector(data_g_t_intervals), 
                     xy_intvl = as.vector(data_g_xy_intervals), 
                     p_ij = as.vector(p_mat))
  g_df <- g_df[g_df$t_intvl > 0 & g_df$xy_intvl > 0]
  g_t = g_df[, list(g_t_fun = sum(p_ij)/offspring_sum), by = t_intvl]
  setkey(g_t, t_intvl)
  g_t$g_t_fun = g_t$g_t_fun/g_t_bw
  g_xy <- g_df[, list(g_xy_fun = sum(p_ij)/offspring_sum), by = xy_intvl]
  setkey(g_xy, xy_intvl)
  g_xy$g_xy_fun = g_xy$g_xy_fun/g_xy_bw
  
  lambda_mat = matrix(0, N, N)
  lam_index <- which(data_g_t_intervals > 0 & data_g_xy_intervals > 0)
  lambda_mat[lam_index] = alpha*g_t$g_t_fun[data_g_t_intervals[lam_index]]*(g_xy$g_xy_fun[data_g_xy_intervals[lam_index]])/(2*pi*g_xy_dist_mat[lam_index])
  diag(lambda_mat) <- mu_vec
  
  p_mat_upd = lambda_mat/rowSums(lambda_mat)
  
  delta = max(abs(p_mat - p_mat_upd))
  print(delta)
  
  p_mat = p_mat_upd
  
  counter = counter + 1
}

p1 = qplot(x = g_t_breaks[-1], y = g_t$g_t_fun, main = 'G_T_FUN', xlab =  'Time', ylab = 'G_T')
p2 = qplot(y = g_xy$g_xy_fun, main = 'G_XY_FUN', xlab =  'Distance', ylab = 'G_XY')

print(ggarrange(p1, p2))

sum(g_t$g_t_fun*g_t_bw)
sum(g_xy$g_xy_fun*g_xy_bw)
setkeyv(mu_txy, c('binx', 'biny'))

#
g_t <- rbind(g_t, cbind(max(g_t$t_intvl) + 1, 0), use.names = F)
g_xy <- rbind(g_xy, cbind(max(g_xy$xy_intvl) + 1, 0), use.names = F)
g_t_breaks <- c(g_t_breaks, 3000)
g_xy_breaks <- c(g_xy_breaks, 3)

# save data
data <- data[,2:4]
save(mu_txy, g_t, g_xy, alpha, mu_x_breaks, mu_y_breaks, mu_df, g_t_breaks, g_xy_breaks, file = 'estimates.Rdata')