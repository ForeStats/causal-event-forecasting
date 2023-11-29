# helper script
InHomoDiscK <- function(s, events, N, p_hat){
  
  events <- events[events <= N + s]
  
  outersum <- sapply(events[events <= N], FUN = function(z){
    
    length(events[events > z & events <= z + s])
  })
  
  sum(outersum)/(N*(p_hat^2))
}

event_simulator <- function(t, params, events){
  
  #Parameters
  mu_txym   <- params[[1]]
  alpha     <- params[[2]]
  g_t       <- params[[3]]
  g_xy      <- params[[4]]
  
  mu_x_breaks    <- params[[5]]
  mu_y_breaks    <- params[[6]]
  g_t_breaks     <- params[[7]]
  g_xy_breaks    <- params[[8]]
  
  gamma_mat      <- params[[9]]
  
  mu_x_bw  <- params[[10]]
  mu_y_bw  <- params[[11]]
  g_t_bw   <- params[[12]]
  g_xy_bw  <- params[[13]]
  
  # Initializations
  data <- data.frame()
  
  # Background events
  mu_txym$bevents <- rpois(nrow(mu_txym), (mu_txym$mu)*(mu_x_bw*mu_y_bw))
  
  bgdata <- data.frame()
  for(i in 1:nrow(mu_txym)){
    if(mu_txym$bevents[i] > 0){
      newx <- runif(mu_txym$bevents[i], mu_x_breaks[mu_txym$binx[i]], mu_x_breaks[mu_txym$binx[i]+1])
      newy <- runif(mu_txym$bevents[i], mu_y_breaks[mu_txym$biny[i]], mu_y_breaks[mu_txym$biny[i]+1])
      
      bgdata <- cbind(t = t, x = newx, y = newy, m = mu_txym$type[i], Z = 0)
    }
  }
  
  # Triggered events
  # Lambda at time t
  t_intvls <- cut(t-events$t, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
  # Need to multiply by t_bw if width > 1
  lambda_t_vec <- alpha*g_t$g_t_fun[t_intvls]
  trigger_events <- rpois(length(lambda_t_vec), lambda_t_vec)
  
  tgdata <- data.frame()
  for(i in which(trigger_events > 0)){
    r_intvl <- sample(g_xy$xy_intvl, trigger_events[i], prob = g_xy$g_xy_fun*g_xy_bw, replace = TRUE)
    r <- runif(trigger_events[i], g_xy_breaks[r_intvl], g_xy_breaks[r_intvl+1])
    theta <- runif(trigger_events[i], 0, 2*pi)
    newx <- events$x[i] + r*cos(theta)
    newy <- events$y[i] + r*sin(theta)
    newm <- sample(1:ncol(gamma_mat), trigger_events[i], prob = gamma_mat[events$m[i],], replace = TRUE)
    
    tgdata <- cbind(t = t, x = newx, y = newy, m = newm, Z = i)  
  }
  
  data <- rbind(data, bgdata, tgdata)
    
  data <- data[(data$x > min(mu_x_breaks) & data$x < max(mu_x_breaks)) & (data$y > min(mu_y_breaks) & data$y < max(mu_y_breaks)),]
  rownames(data) <- NULL
  
  return(data)
}
