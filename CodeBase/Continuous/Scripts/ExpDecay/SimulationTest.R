# Simulation Script
library(Hmisc)
library(distances)
options(digits=10)
#set.seed(123)

# Truth Setup
T_END <- 1000
mu_t_bw <- T_END
mu_x_breaks <- seq(0, 1, 0.5)
mu_y_breaks <- seq(0, 1, 0.5)
mu_x_bw <- diff(mu_x_breaks); mu_y_bw <- diff(mu_y_breaks);
mu_txym <- expand.grid(bint = 1, binx = 1:2, biny = 1:2, binm = 1:4)
mu_txym$mu <- c(rep(0.5, 4), rep(0.1, 4), rep(0.3, 4), rep(1, 4))

g_t_breaks = c(seq(0, 10, 1), T_END)
g_t_bw = diff(g_t_breaks)
g_t <- data.frame(t_intvl = 1:11, g_t_fun = c(exp(seq(0, -18, -2)),0))
g_t$g_t_fun <- g_t$g_t_fun/sum(g_t$g_t_fun*g_t_bw)

g_xy_breaks = c(seq(0, 0.1, 0.01), 2)
g_xy_bw = diff(g_xy_breaks)
g_xy <- data.frame(xy_intvl = 1:11, g_xy_fun = c(exp(seq(0, -18, -2)),0))
g_xy$g_xy_fun <- g_xy$g_xy_fun/sum(g_xy$g_xy_fun*g_xy_bw)

sum(g_t$g_t_fun*g_t_bw)
sum(g_xy$g_xy_fun*g_xy_bw)

gamma_mat <- matrix(c(0.3, 0.3, 0.3, 0.1,
                      0,   0.7, 0,   0.3,
                      0,   0.5, 0.5, 0,
                      0,   0.7, 0.2, 0.1), 4, 4, byrow = T)

alpha = 0.8

# Background events

bgdata <- data.frame()
for(i in 1:nrow(mu_txym)){
  mu_i <- (mu_txym$mu[i])*(mu_x_bw[mu_txym$binx[i]])*(mu_y_bw[mu_txym$biny[i]])
  n_i  <- as.integer(mu_i*T_END)
  newt <- cumsum(rexp(2*n_i, mu_i))
  newt <- newt[newt < T_END]
  n_i <- length(newt)
    
  if(n_i > 0){
    
    newx <- runif(n_i, mu_x_breaks[mu_txym$binx[i]], mu_x_breaks[mu_txym$binx[i]+1])
    newy <- runif(n_i, mu_y_breaks[mu_txym$biny[i]], mu_y_breaks[mu_txym$biny[i]+1])
    
    bgdata <- rbind(bgdata, cbind(t = newt, x = newx, y = newy, m = mu_txym$binm[i], Z = 0))
  }
}

bgdata <- bgdata[order(bgdata$t),]
Immg <- bgdata

#Clusters
for(i in 1:nrow(Immg)){
  
  CurGen <- Immg[i,]
  
  while(nrow(CurGen) > 0){
    
    for(j in 1:nrow(CurGen)){
      
      n_child <- rpois(1, lambda = alpha)
      
      if(n_child > 0){
        
        # Triggered events
        t_intvls <- sample(g_t$t_intvl, n_child, prob = g_t$g_t_fun*g_t_bw, replace = TRUE)
        newt <- CurGen$t[1] + runif(n_child, g_t_breaks[t_intvls], g_t_breaks[t_intvls+1])
        
        r_intvl <- sample(g_xy$xy_intvl, n_child, prob = g_xy$g_xy_fun*g_xy_bw, replace = TRUE)
        r <- runif(n_child, g_xy_breaks[r_intvl], g_xy_breaks[r_intvl+1])
        theta <- runif(n_child, 0, 2*pi)
        newx <- CurGen$x[1] + r*cos(theta)
        newy <- CurGen$y[1] + r*sin(theta)
        
        newm <- sample(1:ncol(gamma_mat), n_child, prob = gamma_mat[CurGen$m[1],], replace = TRUE)
        
        new_events <- data.frame(t = newt, x = newx, y = newy, m = newm, Z = CurGen$t[1])
        bgdata <- rbind(bgdata, new_events)
        CurGen <- rbind(CurGen, new_events)
      }
      CurGen <- CurGen[-1, ]
    }
    CurGen <- CurGen[order(CurGen$t), ]
  }
}

bgdata <- bgdata[order(bgdata$t),]
data <- bgdata[(bgdata$x > 0 & bgdata$x < 1) & (bgdata$y > 0 & bgdata$y < 1) & (bgdata$t < T_END),]
rownames(data) <- NULL

data$B <- vapply(data$Z, FUN = function(z){
  if(z > 0){
    ifelse(length(which(data$t == z)) > 0, which(data$t == z), 0)
  }else{
    0
  }}, FUN.VALUE = 0)
table(data$B > 0)/nrow(data)

train_data <- data[data$t > 100 & data$t < 500, ]
train_data$t <- train_data$t - 100

save(data, train_data, file = 'simdata.RData')
