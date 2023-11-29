# Load packages
library(Hmisc)
library(distances)
library(data.table)
library(ggplot2)
library(ggpubr)

# Load data
setwd("~/Documents/MyResearch/CodeBase/Continuous/Scripts/BellDecay")
load("~/Documents/MyResearch/CodeBase/Continuous/Scripts/BellDecay/simdata.RData")

data <- train_data
rm(train_data)
# times and x,y
T_END <- 400
N <- nrow(data)
data <- cbind(index = 1:N, data)
row.names(data) <- NULL
type_counts <- as.numeric(table(data$m))

# Initialise the P matrix
p_mat <- matrix(1, N, N)
p_mat[upper.tri(p_mat)] <- 0
p_mat <- p_mat/(1:N)

# Algorithm
# Background bandwidth
mu_t_bw <- T_END
mu_x_breaks <- seq(0, 1, 0.5)
mu_y_breaks <- seq(0, 1, 0.5)
mu_x_bw <- diff(mu_x_breaks); mu_y_bw <- diff(mu_y_breaks);

# Set bins
data$bint <- 1
data$binx <- as.numeric(cut2(data$x, cuts = mu_x_breaks, digits = 7))
data$biny <- as.numeric(cut2(data$y, cuts = mu_y_breaks, digits = 7))

# Excitation bandwidth
g_t_dist_mat <- as.matrix(distances(data$t))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

g_t_breaks = seq(0, 10, 1)
g_t_bw = diff(g_t_breaks)

data_g_t_intervals = apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0

g_xy_dist_mat <- as.matrix(distances(data[,c('x', 'y')]))
g_xy_dist_mat[upper.tri(g_xy_dist_mat)] <- NA

g_xy_breaks = seq(0, 0.1, 0.01)
g_xy_bw = diff(g_xy_breaks)

data_g_xy_intervals = apply(g_xy_dist_mat, 2, FUN = cut, breaks = g_xy_breaks, labels = F, include.lowest = TRUE)
data_g_xy_intervals[data_g_t_intervals == 0] <- 0
data_g_xy_intervals[is.na(data_g_xy_intervals)] <- 0

delta = 1
counter = 1
mu_data <- data.table(data[,c('index', 'bint', 'binx', 'biny', 'm')], p_ii = diag(p_mat))
event_types <- as.numeric(data$m)

while(delta > 1e-3){
  # Update mu
  pii_vec <- diag(p_mat)
  mu_data$p_ii <- pii_vec
  mu_txym = mu_data[, list(mu = sum(p_ii)), by = c('bint', 'binx', 'biny', 'm')]
  mu_txym$mu = mu_txym$mu/(mu_t_bw[mu_txym$bint])
  mu_txym$mu = mu_txym$mu/(mu_x_bw[mu_txym$binx])
  mu_txym$mu = mu_txym$mu/(mu_y_bw[mu_txym$biny])
  mu_data = mu_txym[mu_data, on = c('bint', 'binx', 'biny', 'm')]
  setkey(mu_data, index)
  
  # Update triggering function
  offspring_sum <- N - sum(pii_vec)
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
  
  # Gamma mat
  gamma_p_mat <- p_mat
  diag(gamma_p_mat) <- 0
  gamma_p_mat <- cbind(child = event_types, gamma_p_mat)
  colnames(gamma_p_mat)[-1] <- event_types
  gamma_melt <- melt(data.table(gamma_p_mat), id.vars = 'child', variable.name = 'parent', value.name = 'p_ij', variable.factor = FALSE)
  gamma_df = gamma_melt[, list(gamma_fun = sum(p_ij)), by = c('child', 'parent')]
  setorder(gamma_df, parent, child)
  gamma_cast <- as.matrix(dcast.data.table(gamma_df, formula = child ~ parent, value.var = 'gamma_fun'))[,-1]
  gamma_cast <- t(t(gamma_cast)/colSums(gamma_cast))
  gamma_mat <- gamma_cast[event_types , ]
  gamma_mat <- gamma_mat[ , event_types]
  
  lambda_mat = matrix(0, N, N)
  lam_index <- which(data_g_t_intervals > 0 & data_g_xy_intervals > 0)
  lambda_mat[lam_index] = alpha*gamma_mat[lam_index]*g_t$g_t_fun[data_g_t_intervals[lam_index]]*(g_xy$g_xy_fun[data_g_xy_intervals[lam_index]])/(2*pi*g_xy_dist_mat[lam_index])
  diag(lambda_mat) <- mu_data$mu
  
  p_mat_upd = lambda_mat/rowSums(lambda_mat)
  
  delta = max(abs(p_mat - p_mat_upd))
  print(alpha)
  
  p_mat = p_mat_upd
  
  counter = counter + 1
}

p1 = qplot(y = g_t$g_t_fun, main = 'G_T_FUN', xlab =  'Time', ylab = 'G_T')
p2 = qplot(y = g_xy$g_xy_fun, main = 'G_XY_FUN', xlab =  'Distance', ylab = 'G_XY')

print(ggarrange(p1, p2))

round(t(gamma_cast) ,2)
sum(g_t$g_t_fun*g_t_bw)
sum(g_xy$g_xy_fun*g_xy_bw)
setkeyv(mu_txym, c('binx', 'biny', 'm'))
data$predZ <- apply(p_mat, 1, which.max)
data$predZ[data$predZ == data$index] <- 0
table(data$Z > 0)
table(data$predZ > 0)
table(data$Z == 0, data$predZ == 0)
