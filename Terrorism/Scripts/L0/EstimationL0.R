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

T_END <- 7*52*8
data  <- data[time <= T_END & bint == 2]

M     <- nrow(type_meta_data)
N     <- nrow(data)
data  <- cbind(index = 1:N, data)
row.names(data) <- NULL
type_counts     <- as.numeric(table(data$type))

# Background bandwidth
mu_t_bw      <- 7*52*2

# Excitation bandwidths
g_t_dist_mat <- as.matrix(distances(data$time))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

#quantile(g_t_dist_mat[g_t_dist_mat > 0], seq(0, 0.1, 0.01), na.rm = TRUE)
g_t_breaks <- c(0, 1, 2^(3:6)-1)
g_t_bw     <- diff(g_t_breaks)

data_g_t_intervals <- apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0
#table(data_g_t_intervals[data_g_t_intervals>0])

# Initialisations
delta   <- 1
counter <- 1

# Initialise the P matrix
p_mat <- matrix(1, N, N)
p_mat[upper.tri(p_mat)] <- 0
p_mat <- p_mat/(1:N)

# Lambda mat
lambda_mat <- matrix(0, N, N)
lam_index  <- which(data_g_t_intervals > 0)
diag_index <- which(row(lambda_mat) == col(lambda_mat))
lower_tri_index <- sort(c(lam_index, diag_index))

data$bint    <- 1
mu_data      <- data.table(data[,c('index', 'bint', 'type')], p_ii = p_mat[diag_index])
event_types  <- as.numeric(data$type)
gamma_mat    <- alpha_mat <- matrix(event_types, N, N, byrow = T)

# Meta df
meta_p_mat               <- p_mat
meta_p_mat[diag_index]   <- 0
meta_p_mat               <- cbind(child_index = 1:N, meta_p_mat)
colnames(meta_p_mat)[-1] <- 1:N
meta_melt                <- melt(data.table(meta_p_mat), id.vars = 'child_index', variable.name = 'parent_index', 
                                 value.name = 'p_ij', variable.factor = FALSE)
meta_melt$parent_index   <- as.numeric(meta_melt$parent_index)
meta_melt$child_type     <- data$type[meta_melt$child_index]
meta_melt$parent_type    <- data$type[meta_melt$parent_index]
meta_melt$t_intvl        <- as.vector(data_g_t_intervals)
meta_melt                <- meta_melt[lam_index,]

while(delta > 1e-3){
  # Update mu
  mu_data$p_ii <- p_mat[diag_index]
  mu_tm        <- mu_data[, list(mu = sum(p_ii)), by = c('bint', 'type')]
  mu_tm$mu     <- mu_tm$mu/mu_t_bw
  mu_data      <- mu_tm[mu_data, on = c('bint', 'type')]
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
  
  # Gamma function
  gamma_df   <- meta_melt[, list(gamma_fun = sum(p_ij)), by = c('child_type', 'parent_type')]
  setorder(gamma_df, parent_type, child_type)
  gamma_cast <- as.matrix(dcast.data.table(gamma_df, formula = child_type ~ parent_type, value.var = 'gamma_fun'))[,-1]
  gamma_cast <- t(t(gamma_cast)/colsums(gamma_cast))
  meta_melt$gamma <- gamma_cast[cbind(meta_melt$child_type, meta_melt$parent_type)]
  
  # Alpha function
  meta_melt$alpha <- offspring_sum/N
  
  # Lambda Calc
  meta_melt  <- meta_melt[ , lambda := (alpha*gamma*g_t) ]
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

plot(g_t$g_t_fun)

round(t(gamma_cast) ,2)
colSums(g_t[,-1]*g_t_bw)
round(g_t[,-1]*g_t_bw, 2)
setkeyv(mu_tm, c('bint', 'type'))

gamma_plot <- t(gamma_cast)
type_meta_data$short_label <- c('ARMED ASSAULT', 'ASSASSINATION', 'BOMBING', 'FACILITY ATTACK', 'KIDNAPPING')
row.names(gamma_plot)      <- colnames(gamma_plot) <- type_meta_data$short_label

pdf("~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/L0/CorrPlotL0.pdf", width=8, height=8)
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
g_t_breaks  <- c(g_t_breaks, 2^12)
alpha       <- meta_melt$alpha[1]
# save data
save(mu_tm, g_t, alpha, gamma_cast, g_t_breaks,  
     file = '~/Documents/MyResearch/Projects/Terrorism/BiharChat/Scripts/L0/estimatesL0.Rdata')

library(data.tree)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
source('helpers.R')
set.seed(12345)

# Branching structure
data$Z   <- apply(p_mat, 1, FUN = function(z) sample(1:N, 1, prob = z))
data$Z[data$Z == data$index] <- 0
data$Gen <- sapply(data$index, FUN = function(z) findGen(ind = z, data = data))
table(data$Z > 0)
head(sort(table(data$Z), decreasing = T), 10)
tree <- tree_plot(indx = 480, depth = 4, data = data, type_meta_data = type_meta_data)
export_graph(ToDiagrammeRGraph(tree), "~/Documents/MyResearch/Projects/Terrorism/BiharChat/Figures/L0/treeL0.pdf")
