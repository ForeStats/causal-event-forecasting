# Simulation Script
library(Hmisc)
library(distances)
options(digits=10)
set.seed(123)

# Truth Setup
T_END <- 10000
mu_t_bw <- T_END
mu_x_bw <- 50; mu_y_bw  <- 50
mu_x_breaks <- seq(0, 100, mu_x_bw)
mu_y_breaks <- seq(0, 50, mu_y_bw)
mu_txym <- expand.grid(bint = 1, binx = 1:2, biny = 1, binm = 1:2)
mu_txym$mu <- c(0.05, 0.01, 0.05, 0.01)/(mu_x_bw*mu_y_bw)

g_t_breaks = c(seq(0, 10, 1), T_END)
g_t_bw = diff(g_t_breaks)
g_t <- data.frame(t_intvl = 1:11, g_t_fun = c(exp(seq(0, -18, -2)),0))
g_t$g_t_fun <- g_t$g_t_fun/sum(g_t$g_t_fun*g_t_bw)

g_xy_breaks = c(seq(0, 10, 1), 200)
g_xy_bw = diff(g_xy_breaks)
g_xy <- data.frame(xy_intvl = 1:11, g_xy_fun = c(exp(seq(0, -18, -2)),0))
g_xy$g_xy_fun <- g_xy$g_xy_fun/sum(g_xy$g_xy_fun*g_xy_bw)

sum(g_t$g_t_fun*g_t_bw)
sum(g_xy$g_xy_fun*g_xy_bw)

gamma_mat <- matrix(c(0.7, 0.3, 0, 1), 2, 2, byrow = T)

alpha = 0.8

# Initializations
data <- data.frame()
t <- 1

while(t < T_END){
  
  # Background events
  mu_txym$bevents <- rpois(nrow(mu_txym), (mu_txym$mu)*(mu_x_bw*mu_y_bw))
  
  bgdata <- data.frame()
  for(i in 1:nrow(mu_txym)){
    if(mu_txym$bevents[i] > 0){
      newx <- runif(mu_txym$bevents[i], mu_x_breaks[mu_txym$binx[i]], mu_x_breaks[mu_txym$binx[i]+1])
      newy <- runif(mu_txym$bevents[i], mu_y_breaks[mu_txym$biny[i]], mu_y_breaks[mu_txym$biny[i]+1])
      
      bgdata <- cbind(t = t, x = newx, y = newy, m = mu_txym$binm[i], Z = 0)
    }
  }
  
  # Triggered events
  # Lambda at time t
  t_intvls <- cut(t-data$t, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
  # Need to multiply by t_bw if width > 1
  lambda_t_vec <- alpha*g_t$g_t_fun[t_intvls]
  trigger_events <- rpois(length(lambda_t_vec), lambda_t_vec)
  
  tgdata <- data.frame()
  for(i in which(trigger_events > 0)){
    r_intvl <- sample(g_xy$xy_intvl, trigger_events[i], prob = g_xy$g_xy_fun*g_xy_bw, replace = TRUE)
    r <- runif(trigger_events[i], g_xy_breaks[r_intvl], g_xy_breaks[r_intvl+1])
    theta <- runif(trigger_events[i], 0, 2*pi)
    newx <- data$x[i] + r*cos(theta)
    newy <- data$y[i] + r*sin(theta)
    newm <- sample(1:ncol(gamma_mat), trigger_events[i], prob = gamma_mat[data$m[i],], replace = TRUE)
    
    tgdata <- cbind(t = t, x = newx, y = newy, m = newm, Z = i)  
  }
  
  data <- rbind(data, bgdata, tgdata)
  
  t <- t+1
}
data <- data[(data$x > 0 & data$x < 100) & (data$y > 0 & data$y < 50),]
rownames(data) <- NULL

name.check <- function(vec, k){
  
  nam <- names(vec)
  add.k   <- which(!1:k %in% as.numeric(names(vec)))
  vec <- append(vec,rep(0,length(add.k)))
  names(vec) <- c(nam,add.k)
  vec <- vec[order(as.numeric(names(vec)))]
  names(vec) <- paste0('T',names(vec))
  
  return(vec)
}
mean(name.check(table(data$Z[data$Z > 0]), nrow(data)))

# PMat
data$bint <- 1
data$binx <- as.numeric(cut2(data$x, cuts = mu_x_breaks, digits = 7))
data$biny <- as.numeric(cut2(data$y, cuts = mu_y_breaks, digits = 7))
table(data$binx, data$biny)

g_t_dist_mat <- as.matrix(distances(data$t))
g_t_dist_mat[upper.tri(g_t_dist_mat)] <- NA

data_g_t_intervals = apply(g_t_dist_mat, 2, FUN = cut, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
data_g_t_intervals[is.na(data_g_t_intervals)] <- 0
table(data_g_t_intervals)

g_xy_dist_mat <- as.matrix(distances(data[,c('x', 'y')]))
g_xy_dist_mat[upper.tri(g_xy_dist_mat)] <- NA

data_g_xy_intervals = apply(g_xy_dist_mat, 2, FUN = cut, breaks = g_xy_breaks, labels = F, include.lowest = TRUE)
data_g_xy_intervals[data_g_t_intervals == 0] <- 0
data_g_xy_intervals[is.na(data_g_xy_intervals)] <- 0
table(data_g_xy_intervals)

gamma_values <- as.matrix(t(gamma_mat))[data$m , ]
gamma_values <- gamma_values[ , data$m]
gamma_values[upper.tri(gamma_values)] <- 0
diag(gamma_values) <- 0
table(gamma_values)

# Initialise the P matrix
N <- nrow(data)
p_mat <- matrix(0, N, N)
p_mat[1,1] <- 1

for(i in 2:N){
  # Background
  mu <- mu_txym$mu[mu_txym$binx == data$binx[i] & mu_txym$biny == data$biny[i] & mu_txym$binm == data$m[i]]
  lam_index <- which(data_g_t_intervals[i,] > 0 & data_g_xy_intervals[i,] > 0)
  lambda_vec = (alpha*gamma_values[i,][lam_index]*g_t$g_t_fun[data_g_t_intervals[i,][lam_index]]*(g_xy$g_xy_fun[data_g_xy_intervals[i,][lam_index]]))/(2*pi*(g_xy_dist_mat[i,][lam_index]))
  lambda <- c(lambda_vec, mu)
  p_mat[i, 1:length(lambda)] <- lambda/sum(lambda)
}
sum(diag(p_mat))
table(data$Z > 0)
gamma_p_mat <- p_mat
diag(gamma_p_mat) <- 0
mean(colSums(gamma_p_mat))
(N - sum(diag(p_mat)))/N

train_data <- data[data$t > 1000 & data$t < 5000, ]
train_data$t <- train_data$t - 1000

save(data, train_data, file = 'simdata.RData')
