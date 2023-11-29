# Load packages
library(Hmisc)
library(distances)
library(data.table)
library(Rfast)
library(ggplot2)
library(ggpubr)
library(corrplot)

# Load data
setwd("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts")
load("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/data.RData")

#T_END <- 56*52
#data  <- data[time <= T_END]

T_END <- 2000
data  <- data[time > (max(data$time) - 2000)]

M     <- nrow(type_meta_data)
N     <- nrow(data)
data  <- cbind(index = 1:N, data)
row.names(data) <- NULL
type_counts     <- as.numeric(table(data$type))

# Background bandwidth
#mu_t_bw     <- T_END
mu_t_bw     <- T_END
mu_x_breaks <- seq(longlims[1], longlims[2], diff(longlims)/17)
mu_y_breaks <- seq(latlims[1], latlims[2], diff(latlims)/17)
mu_x_bw     <- diff(mu_x_breaks); mu_y_bw <- diff(mu_y_breaks);

# Set bins
data$bint <- 1
data$binx <- as.numeric(cut2(data$longitude, cuts = mu_x_breaks, digits = 7))
data$biny <- as.numeric(cut2(data$latitude, cuts = mu_y_breaks, digits = 7))

# Excitation bandwidths
g_t_dist_mat <- as.matrix(distances(data$time))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

#quantile(g_t_dist_mat[g_t_dist_mat > 0], seq(0, 0.1, 0.01), na.rm = TRUE)
g_t_breaks <- c(0, 1, (2^(2:6)-1))
g_t_bw     <- diff(g_t_breaks)

data_g_t_intervals <- apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0
#table(data_g_t_intervals[data_g_t_intervals>0])

g_xy_dist_mat <- as.matrix(distances(data[,c('longitude', 'latitude')]))
g_xy_dist_mat[upper.tri(g_xy_dist_mat)] <- NA
#min(g_xy_dist_mat[g_xy_dist_mat > 0], na.rm = TRUE)
g_xy_dist_mat[g_xy_dist_mat == 0 & is.finite(g_xy_dist_mat)] <- 5e-7

g_xy_breaks <- c(0, 1e-6, 0.1, 0.2)
g_xy_bw     <- diff(g_xy_breaks)

data_g_xy_intervals <- apply(g_xy_dist_mat, 2, FUN = cut, breaks = g_xy_breaks, labels = F, include.lowest = TRUE)
data_g_xy_intervals[data_g_t_intervals == 0] <- 0
data_g_xy_intervals[is.na(data_g_xy_intervals)] <- 0
#table(data_g_xy_intervals[data_g_xy_intervals>0])

# Initialisations
delta   <- 1
counter <- 1

# Initialise the P matrix
p_mat <- matrix(1, N, N)
p_mat[upper.tri(p_mat)] <- 0
p_mat <- p_mat/(1:N)

# Lambda mat
lambda_mat <- matrix(0, N, N)
lam_index  <- which(data_g_t_intervals > 0 & data_g_xy_intervals > 0)
diag_index <- which(row(lambda_mat) == col(lambda_mat))
lower_tri_index <- sort(c(lam_index, diag_index))

mu_data <- data.table(data[,c('index', 'bint', 'binx', 'biny', 'type')], p_ii = p_mat[diag_index])
event_types <- as.numeric(data$type)
gamma_mat   <- alpha_mat <- matrix(event_types, N, N, byrow = T)

# Meta df
meta_p_mat       <- p_mat
meta_p_mat[diag_index]   <- 0
meta_p_mat       <- cbind(child_index = 1:N, meta_p_mat)
colnames(meta_p_mat)[-1] <- 1:N
meta_melt                <- melt(data.table(meta_p_mat), id.vars = 'child_index', variable.name = 'parent_index', 
                                 value.name = 'p_ij', variable.factor = FALSE)
meta_melt$parent_index   <- as.numeric(meta_melt$parent_index)
meta_melt$child_type     <- data$type[meta_melt$child_index]
meta_melt$parent_type    <- data$type[meta_melt$parent_index]
meta_melt$t_intvl        <- as.vector(data_g_t_intervals)
meta_melt$xy_intvl       <- as.vector(data_g_xy_intervals)
meta_melt                <- meta_melt[lam_index,]
meta_melt$xy_dist        <- 2*pi*g_xy_dist_mat[lam_index]

while(delta > 1e-4){
  # Update mu
  mu_data$p_ii <- p_mat[diag_index]
  mu_txym      <- mu_data[, list(mu = sum(p_ii)), by = c('bint', 'binx', 'biny', 'type')]
  mu_txym$mu   <- mu_txym$mu/((mu_t_bw[mu_txym$bint])*(mu_x_bw[mu_txym$binx])*(mu_y_bw[mu_txym$biny]))
  mu_data      <- mu_txym[mu_data, on = c('bint', 'binx', 'biny', 'type')]
  setkey(mu_data, index)
  
  # Update triggering functions
  meta_melt$p_ij <- p_mat[lam_index]
  
  # Offspring function
  offspring_sum <- N - sum(mu_data$p_ii)
  
  # G(T) function
  g_t <- meta_melt[, list(g_t_fun = sum(p_ij)/offspring_sum), by = t_intvl]
  setorder(g_t, t_intvl)
  g_t$g_t_fun   <- g_t$g_t_fun/g_t_bw
  meta_melt$g_t <- g_t$g_t_fun[meta_melt$t_intvl]
  
  # H(XY) function
  g_xy <- meta_melt[, list(g_xy_fun = sum(p_ij)/offspring_sum), by = xy_intvl]
  setorder(g_xy, xy_intvl)
  g_xy$g_xy_fun   <- g_xy$g_xy_fun/g_xy_bw
  meta_melt$g_xy  <- g_xy$g_xy_fun[meta_melt$xy_intvl]
  
  # Gamma function
  gamma_df   <- meta_melt[, list(gamma_fun = sum(p_ij)), by = c('child_type', 'parent_type')]
  setorder(gamma_df, parent_type, child_type)
  gamma_cast <- as.matrix(dcast.data.table(gamma_df, formula = child_type ~ parent_type, value.var = 'gamma_fun'))[,-1]
  gamma_cast <- t(t(gamma_cast)/colsums(gamma_cast))
  meta_melt$gamma <- gamma_cast[cbind(meta_melt$child_type, meta_melt$parent_type)]
  
  # Alpha function
  meta_melt$alpha <- offspring_sum/N
  
  # Lambda Calc
  meta_melt  <- meta_melt[ , lambda := (alpha*gamma*g_t*g_xy)/(xy_dist) ]
  lambda_mat[lam_index]  <- meta_melt$lambda
  lambda_mat[diag_index] <- mu_data$mu
  mu_data$mu <- NULL
  
  # Update Pij
  p_mat_upd <- lambda_mat/rowsums(lambda_mat)
  delta     <- max(abs( p_mat[lower_tri_index] - p_mat_upd[lower_tri_index] ))
  p_mat     <- p_mat_upd
  
  # Print progress
  print(paste0("Iter: ", counter, " alpha: ", round(offspring_sum/N, 4), " delta: ", round(delta, 6)))
  counter   <- counter + 1
}

p1 <- qplot(x = log(g_t_breaks[-1]), y = log(g_t$g_t_fun), main = 'G_T_FUN', xlab =  'Time', ylab = 'G_T')
p2 <- qplot(x = log(g_xy_breaks[-1]), y = log(g_xy$g_xy_fun), main = 'G_XY_FUN', xlab =  'Distance', ylab = 'G_XY')

print(ggarrange(p1, p2))
plot(g_t$g_t_fun)
plot(g_xy$g_xy_fun)

round(t(gamma_cast) ,2)
colSums(g_t[,-1]*g_t_bw)
round(g_t[,-1]*g_t_bw, 2)
colSums(g_xy[,-1]*g_xy_bw)
round(g_xy[,-1]*g_xy_bw, 2)
setkeyv(mu_txym, c('binx', 'biny', 'type'))

gamma_plot <- t(gamma_cast)
type_meta_data$short_label <- c('ARMED ASSAULT', 'ASSASSINATION', 'BOMBING', 'FACILITY ATTACK', 'KIDNAPPING')
row.names(gamma_plot)      <- colnames(gamma_plot) <- type_meta_data$short_label

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/CorrPlot.pdf", width=8, height=8)
corrplot(gamma_plot,
         method="circle",
         is.corr=FALSE,
         type="full",
         tl.srt = 45,
         number.cex = 1,
         cl.lim=c(0,1), 
         addCoef.col = rgb(0,0,0, alpha = 0.6)
)
dev.off()

# Make sure extends validation time as well
g_t  <- rbind(g_t, cbind(max(g_t$t_intvl) + 1, 0), use.names = F)
g_xy <- rbind(g_xy, cbind(max(g_xy$xy_intvl) + 1, 0), use.names = F)
g_t_breaks  <- c(g_t_breaks, 2^12)
g_xy_breaks <- c(g_xy_breaks, 17)
alpha       <- meta_melt$alpha[1]
# save data
save(mu_txym, g_t, g_xy, alpha, gamma_cast, mu_x_breaks, mu_y_breaks, g_t_breaks, g_xy_breaks, 
     file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/estimates.Rdata')

library(data.tree)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
source('helpers.R')
# Branching structure
# dim(p_mat)
data$Z <- apply(p_mat, 1, FUN = function(z) which.max(z))
data$Z[data$Z == data$index] <- 0
data$Gen <- sapply(data$index, FUN = function(z) findGen(ind = z, data = data))

head(sort(table(data$Z), decreasing = T),10)
tree <- tree_plot(indx = 563, depth = 10, data = data, type_meta_data = type_meta_data)
export_graph(ToDiagrammeRGraph(tree), "~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/tree.pdf")
