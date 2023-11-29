# helper script
InHomoDiscK <- function(s, events, N, p_hat){
  
  events <- events[events <= N + s]
  
  outersum <- vapply(events[events <= N], FUN = function(z){
    
    length(events[events > z & events <= z + s])
  }, FUN.VALUE = 10)
  
  sum(outersum)/(N*(p_hat^2))
}

findGen <- function(indx, data) {
  if (data$Z[indx] == 0) {
    return(0)
  } else{
    return(1 + findGen(data$Z[indx], data))
  }
}

next.leaf <- function(node){
  sib <- node$siblings
  #print(node$parent$name)
  temp <- node$parent
  if(length(sib)>0){
    for (i in c(1:length(sib))){
      node.sib <- sib[[i]]
      if(node.sib$isLeaf){
        temp <- node.sib
      }
    }
  }
  return(temp)
}

tree_plot <- function(indx, depth, data, type_meta_data){
  
  parent <- Node$new(paste0(data$date[indx], ", ", type_meta_data$short_label[data$type[indx]]))
  parent$index <- indx
  children <- which(data$Z == parent$index)
  if (length(children > 0)){
    for (i in children){
      name.new <- paste0(data$date[i], ", ", type_meta_data$short_label[data$type[i]])
      parent$AddChild(name.new, index = i)
    }
  }else{ 
    parent$AddChild("X")
  }
  
  node <- parent$children[[1]]
  if(node$name == "X"){
    node <- parent
  }
  
  while (node$name != parent$name){
    #print("i")
    if (node$level>depth||node$name == "X"){
      #print("j")
      node <- node$parent
      while(!(node$isLeaf) & node$name != parent$name){
        node <- next.leaf(node)
      }
    }
    if (node$level<=depth & node$isLeaf & node$name != parent$name){
      children <- which(data$Z == node$index)
      if (length(children > 0)){
        for (i in children){
          name.new <- paste0(data$date[i], ", ", type_meta_data$short_label[data$type[i]])
          node$AddChild(name.new, index = i)
        }
      }else{ 
        node$AddChild("X")
      }
      if (node$leafCount > 0){
        node <- node$children[[1]]
      }
    }
  }

  Prune(parent, function(x) x$name != 'X')
  SetGraphStyle(parent, rankdir = "TB")
  SetEdgeStyle(parent, arrowhead = "vee", color = "grey35", penwidth = 2)
  SetNodeStyle(parent, style = "filled,rounded", shape = "box", fillcolor = "LightBlue", 
               fontname = "helvetica", fontcolor = "Black", tooltip = GetDefaultTooltip)
  
  return(parent)
}

event_simulator_L0 <- function(start_t, t_width, params, events){
  
  #Parameters
  mu_tm       <- params[[1]]
  alpha       <- params[[2]]
  g_t         <- params[[3]]
  g_t_breaks  <- params[[4]]
  gamma_mat   <- params[[5]]
  g_t_bw      <- params[[6]]
  
  # Initializations
  events <- events[events$t > (start_t + t_width - tail(g_t_breaks, 2)[1]), ]
  
  data <- data.frame()
  for(j in 1:t_width){
    
    bgdatarows <- apply(mu_tm, 1, FUN = function(z) {
      n_i  <- rpois(1, (z[['mu']]) )
      
      if(n_i > 0){
        return(cbind(t = start_t + j, m = z[['type']]))
      }
    })
    bgdata <- data.frame()
    if(!is.null(bgdatarows)){
      bgdata = Reduce(rbind, bgdatarows)
    }
    
    # Triggered events
    # Lambda at time t
    t_intvls <- cut(start_t + j - events$t, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
    # Need to multiply by t_bw if width > 1
    lambda_t_vec <- alpha*g_t$g_t_fun[t_intvls]
    trigger_events <- rpois(length(lambda_t_vec), lambda_t_vec)
    
    tgdata <- data.frame()
    if(length(which(trigger_events > 0)) > 0){
      tgdata = lapply(which(trigger_events > 0), FUN = function(z) {
        newm <- sample(1:ncol(gamma_mat), trigger_events[z], prob = gamma_mat[events$m[z],], replace = TRUE)
        
        return(cbind(t = start_t + j, m = newm))
      })
      tgdata = Reduce(rbind, tgdata)
    }
    
    if(nrow(bgdata) > 0 & ncol(bgdata) == 2){
      events <- rbind(events, bgdata)
      data   <- rbind(data, bgdata)
    }
    if(nrow(tgdata) > 0 & ncol(tgdata) == 2){
      events <- rbind(events, tgdata)
      data   <- rbind(data, tgdata)
    }
  }
  rownames(data) <- NULL
  
  return(data)
}

event_simulator_L1 <- function(start_t, t_width, params, events){
  
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
  events <- events[events$t > (start_t + t_width - tail(g_t_breaks, 2)[1]), ]
  
  data <- data.frame()
  
  for(j in 1:t_width){
    
    bgdatarows <- apply(mu_txym, 1, FUN = function(z) {
      n_i  <- rpois(1, (z[['mu']])*(mu_x_bw[z[['binx']]])*(mu_y_bw[z[['biny']]]) )
      
      if(n_i > 0){
        newx <- runif(n_i, mu_x_breaks[z[['binx']]], mu_x_breaks[z[['binx']]+1])
        newy <- runif(n_i, mu_y_breaks[z[['biny']]], mu_y_breaks[z[['biny']]+1])
        
        return(cbind(t = start_t + j, x = newx, y = newy, m = z[['type']]))
      }
    })
    bgdata <- data.frame()
    if(!is.null(bgdatarows)){
      bgdata = Reduce(rbind, bgdatarows)
    }
    
    # Triggered events
    # Lambda at time t
    t_intvls <- cut(start_t + j - events$t, breaks = g_t_breaks, labels = F, include.lowest = FALSE)
    # Need to multiply by t_bw if width > 1
    lambda_t_vec <- alpha*g_t$g_t_fun[t_intvls]
    trigger_events <- rpois(length(lambda_t_vec), lambda_t_vec)
    length(which(trigger_events > 0))
    
    tgdata <- data.frame()
    if(length(which(trigger_events > 0)) > 0){
      tgdata = lapply(which(trigger_events > 0), FUN = function(z) {
        r_intvl <- sample(g_xy$xy_intvl, trigger_events[z], prob = g_xy$g_xy_fun*g_xy_bw, replace = TRUE)
        r <- runif(trigger_events[z], g_xy_breaks[r_intvl], g_xy_breaks[r_intvl+1])
        theta <- runif(trigger_events[z], 0, 2*pi)
        newx <- events$x[z] + r*cos(theta)
        newy <- events$y[z] + r*sin(theta)
        newm <- sample(1:ncol(gamma_mat), trigger_events[z], prob = gamma_mat[events$m[z],], replace = TRUE)
        
        return(cbind(t = start_t + j, x = newx, y = newy, m = newm))
      })
      tgdata = Reduce(rbind, tgdata)
    }
    
    data <- rbind(data, bgdata, tgdata)
    if(nrow(bgdata) > 0){
      events <- rbind(events, bgdata)
    }
    if(nrow(tgdata) > 0){
      events <- rbind(events, tgdata)
    }
  }
  
  data <- data[(data$x > min(mu_x_breaks) & data$x < max(mu_x_breaks)) & (data$y > min(mu_y_breaks) & data$y < max(mu_y_breaks)),]
  rownames(data) <- NULL
  
  return(data)
}

calculate_PAI_2d <- function(truth_array, prob_array, zone_prop = NULL, theta = 0, zeta = 0){
  
  hits_spots = rowsums(sapply(1:Z, FUN = function(z){
    rowsums(sapply(1:M, FUN = function(m){
      truth_vec  = truth_array[,z,m]
      
      prob_vec = prob_array[,m]*zone_prop[z]
      
      prob_vec  = 1 - exp(-prob_vec)
      prob_diff = c(0, diff(prob_vec))
      
      hits  = sum(truth_vec[prob_vec > theta | prob_diff > zeta])
      spots = sum(prob_vec > theta | prob_diff > zeta)
      
      return(c(hits, spots))
    }))
  }))
  
  return(hits_spots)
}

calculate_PAI_3d <- function(truth_array, prob_array, zone_prop = NULL, mark_prop = NULL, theta = 0, zeta = 0){
  
  hits_spots = rowsums(sapply(1:Z, FUN = function(z){
    rowsums(sapply(1:M, FUN = function(m){
      truth_vec  = truth_array[,z,m]
      
      prob_vec = prob_array*zone_prop[z]*mark_prop[m]
      
      prob_vec  = 1 - exp(-prob_vec)
      prob_diff = c(0, diff(prob_vec))
      
      hits  = sum(truth_vec[prob_vec > theta | prob_diff > zeta])
      spots = sum(prob_vec > theta | prob_diff > zeta)
      
      return(c(hits, spots))
    }))
  }))
  
  return(hits_spots)
}
