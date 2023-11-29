# Load packages
library(Hmisc)
library(distances)
library(data.table)
library(ggplot2)
library(ggpubr)
library(corrplot)

# Load data
setwd("~/Documents/MyResearch/CodeBase/Discrete/Scripts")
load("~/Documents/MyResearch/CodeBase/Discrete/Workspaces/data.RData")

# times and x,y
T_END <- 2190
data <- data.table(data[data$time <= T_END, ])
N <- nrow(data)
data <- cbind(index = 1:N, data)
row.names(data) <- NULL
type_counts <- as.numeric(table(data$type))

# Initialise the P matrix
p_mat <- matrix(1, N, N)
p_mat[upper.tri(p_mat)] <- 0
p_mat <- p_mat/(1:N)

# Algorithm
# Background bandwidth
range(data$longitude)
range(data$latitude)
mu_t_bw <- T_END
mu_x_bw <- 4; mu_y_bw  <- 4
mu_x_breaks <- seq(83.5, 87.5, mu_x_bw)
mu_y_breaks <- seq(22, 26, mu_y_bw)

data$bint <- 1
data$binx <- as.numeric(cut2(data$longitude, cuts = mu_x_breaks, digits = 7))
data$biny <- as.numeric(cut2(data$latitude, cuts = mu_y_breaks, digits = 7))

# Excitation bandwidth
g_t_dist_mat <- as.matrix(distances(data$time))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

quantile(g_t_dist_mat, seq(0, 1, 0.1), na.rm = T)
#g_t_breaks = seq(0, 100, 1)
#g_t_breaks = c(seq(0, 7, 1), seq(14, 98, 7), seq(126, 500, 28))
g_t_breaks = c(c(0, 1, 2), seq(7, 28, 7), c(56, 500, 1000))
g_t_bw = diff(g_t_breaks)

data_g_t_intervals = apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0

g_xy_dist_mat <- as.matrix(distances(data[,c('longitude', 'latitude')]))
g_xy_dist_mat[upper.tri(g_xy_dist_mat)] <- NA

quantile(g_xy_dist_mat, seq(0, 1, 0.1), na.rm = T)
#g_xy_breaks = seq(0, 1, 0.01)
g_xy_breaks = c(0, rev(exp(-c(0:4, 13:14))), 3, 10)
g_xy_bw = diff(g_xy_breaks)

data_g_xy_intervals = apply(g_xy_dist_mat, 2, FUN = cut, breaks = g_xy_breaks, labels = F, include.lowest = TRUE)
data_g_xy_intervals[data_g_t_intervals == 0] <- 0
data_g_xy_intervals[is.na(data_g_xy_intervals)] <- 0
table(data_g_xy_intervals)

delta = 1
counter = 1
mu_data <- data.table(data[,c('index', 'bint', 'binx', 'biny', 'type')], p_ii = diag(p_mat))
event_types <- as.numeric(data$type)

while(delta > 1e-4){
  # Update mu
  pii_vec <- diag(p_mat)
  mu_data$p_ii <- pii_vec
  mu_txym = mu_data[, list(mu = sum(p_ii)/(mu_t_bw*mu_x_bw*mu_y_bw)), by = c('bint', 'binx', 'biny', 'type')]
  setkeyv(mu_txym, c('bint', 'binx', 'biny', 'type'))
  mu_data = mu_txym[mu_data, on = c('bint', 'binx', 'biny', 'type')]
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
  gamma_cast <- as.matrix(dcast.data.table(gamma_df, formula = child ~ parent, value.var = 'gamma_fun'))[,-1]/type_counts
  gamma_cast <- t(t(gamma_cast)/colSums(gamma_cast))
  gamma_mat <- gamma_cast[event_types , ]
  gamma_mat <- gamma_mat[ , event_types]
  
  lambda_mat = matrix(0, N, N)
  lam_index <- which(data_g_t_intervals > 0 & data_g_xy_intervals > 0)
  lambda_mat[lam_index] = alpha*gamma_mat[lam_index]*g_t$g_t_fun[data_g_t_intervals[lam_index]]*(g_xy$g_xy_fun[data_g_xy_intervals[lam_index]])/(2*pi*g_xy_dist_mat[lam_index])
  diag(lambda_mat) <- mu_data$mu
  
  p_mat_upd = lambda_mat/rowSums(lambda_mat)
  
  delta = max(abs(p_mat - p_mat_upd))
  print(delta)
  
  p_mat = p_mat_upd
  
  counter = counter + 1
}

round(t(gamma_cast) ,2)
sum(g_t$g_t_fun*g_t_bw)
sum(g_xy$g_xy_fun*g_xy_bw)
setkeyv(mu_txym, c('binx', 'biny', 'type'))

plot_data <- data.frame(x = g_t_breaks[-1], y = g_t$g_t_fun)
p1 <- ggplot(plot_data, aes(x=x, y=y)) + 
  geom_point()  +
 # geom_smooth(method = "loess", formula = y ~ x, se=FALSE, col='red', size=0.5, span = 1) +
  xlab("days") + ylab("density") +
  theme_bw()

a <- g_xy_breaks[1:8]
plot_data <- data.frame(x = a[-length(a)] + diff(a)/2, y = g_xy$g_xy_fun[1:7])
p2 <- ggplot(plot_data, aes(x=x, y=y)) + 
  geom_point()  +
 # geom_smooth(method = "loess", formula = y ~ x, se=FALSE, col='red', size=0.5, span = 1) +
  xlab("days") + ylab("density") +
  theme_bw()

print(ggarrange(p1, p2))

gamma_plot <- t(gamma_cast)
type_meta_data <- type_meta_data[order(type_meta_data$type_label),]
row.names(gamma_plot) <- colnames(gamma_plot) <- type_meta_data$type_label

corrplot(gamma_plot,
         method="circle",
         is.corr=FALSE,
         type="full",
         tl.srt = 45,
         number.cex = 1,
         cl.lim=c(0,1), 
         addCoef.col = rgb(0,0,0, alpha = 0.6)
)

#
g_t <- rbind(g_t, cbind(max(g_t$t_intvl) + 1, 0), use.names = F)
g_xy <- rbind(g_xy, cbind(max(g_xy$xy_intvl) + 1, 0), use.names = F)
# Make sure extends validation time as well
g_t_breaks <- c(g_t_breaks, 3000)
g_xy_breaks <- c(g_xy_breaks, 6)

# save data
save(mu_txym, g_t, g_xy, alpha, gamma_cast, mu_x_breaks, mu_y_breaks, g_t_breaks, g_xy_breaks, 
     file = '~/Documents/MyResearch/CodeBase/Discrete/Workspaces/estimates.Rdata')
