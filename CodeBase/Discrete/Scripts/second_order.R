library(Hmisc)
library(MASS)
library(data.table)
library(LaplacesDemon)
library(fitdistrplus)

# Script for second-order analysis
setwd("~/Documents/MyResearch/Terrorism/Code")
load("~/Documents/MyResearch/Terrorism/Code/workspace.RData")

# Projection set-up
R <- 6375
phi0 <- (pi/180)*(sum(range(spp_data$latitude))/2)

spp_data$x <- R*cos(phi0)*(pi/180)*spp_data$longitude
spp_data$y <- R*(pi/180)*spp_data$latitude

# dist of locations
range(spp_data$x)
range(spp_data$y)
spp_data$x <- (spp_data$x - 7120)/100
spp_data$y <- (spp_data$y - 900)/100

# de-duplicates
spp_data <- spp_data[, c('time', 'x', 'y', 'attacktype1')]
spp_data <- unique(spp_data)

# Attack type
# 1 ASSASSINATION
# 2 ARMED ASSAULT
# 3 BOMBING/EXPLOSION
# 4 KIDNAPPING
# 5 FACILITY/INFRASTRUCTURE ATTACK
# 6 OTHERS

spp_data$type_label <- 'OTHERS'
spp_data$type_label[spp_data$attacktype1 == 1] <- 'ASSASSINATION'
spp_data$type_label[spp_data$attacktype1 == 2] <- 'ARMED ASSAULT'
spp_data$type_label[spp_data$attacktype1 == 3] <- 'BOMBING/EXPLOSION'
spp_data$type_label[spp_data$attacktype1 == 6] <- 'KIDNAPPING'
spp_data$type_label[spp_data$attacktype1 == 7] <- 'FACILITY/INFRASTRUCTURE ATTACK'

spp_data$type <- factor(spp_data$type_label, labels = 1:6)

table(spp_data$attacktype1)
table(spp_data$type_label)
table(spp_data$type)

type_meta_data <- unique(spp_data[, c('type', 'type_label')])
spp_data <- spp_data[, c(1,6,2,3)]

# Bins
spp_data$xbin <- Hmisc::cut2(spp_data$x, cuts = seq(from = 0, to = 29, length.out = 30))
table(spp_data$xbin)
spp_data$ybin <- Hmisc::cut2(spp_data$y, cuts = seq(from = 0, to = 30, length.out = 31))
table(spp_data$ybin)
locmatrix <- as.data.frame.matrix(table(spp_data$ybin, spp_data$xbin))
locmatrix <- locmatrix[rev(row.names(locmatrix)),]
#locmatrix[locmatrix < 60] <- 0
#length(locmatrix[locmatrix > 0])

data <- spp_data[, c('time', 'type', 'x', 'y', 'xbin', 'ybin')]
data$diff <- c(1, diff(data$time))
event_times = data$diff
data$diff <- NULL
hist(event_times, right = F, n = 6)

# dist of days to next event 
hist(event_times[event_times > 0], right = F, n = 5)
table(event_times[event_times > 0])/length(event_times[event_times > 0])

xx  <- 0:5
obs <- event_times[event_times > 0] - 1
counts <- as.vector(table(factor(obs, levels = 0:5), exclude = NULL))

poisson.density <- length(obs)*dpois(xx,mean(obs))
nb <- fitdistr(obs,"negative binomial")
nb.density <- length(obs)*dnbinom(xx,size=nb$estimate["size"],mu=nb$estimate["mu"])
geo <- fitdistr(obs,"geometric")
geo.density <- length(obs)*dgeom(xx,prob=geo$estimate["prob"])

foo <- barplot(counts,names.arg=xx,ylim=range(c(counts,poisson.density)))
lines(foo[,1],poisson.density,lwd=2)
lines(foo[,1],nb.density,lwd=2,col="red")
lines(foo[,1],geo.density,lwd=2,col="blue")
legend("topright",lwd=2,col=c("black","red", "blue"),legend=c("Poisson","Negative Binomial", "Geometric"))

# dist of no. of events in day
vec <- rep(0, max(data$time))
names(vec) <- 1:max(data$time)
vec[names(table(data$time))] <- unlist(table(data$time))
hist(vec, n = 25, right = F)
original_distribution = table(vec)/length(vec)

xx  <- 0:13
obs <- vec[vec < 14]
counts <- as.vector(table(factor(obs, levels = 0:13), exclude = NULL))

poisson.density <- length(obs)*dpois(xx,mean(obs))
nb <- fitdistr(obs,"negative binomial")
nb.density <- length(obs)*dnbinom(xx,size=nb$estimate["size"],mu=nb$estimate["mu"])
geo <- fitdistr(obs,"geometric")
geo.density <- length(obs)*dgeom(xx,prob=geo$estimate["prob"])

foo <- barplot(counts,names.arg=xx,ylim=range(c(counts,poisson.density)))
lines(foo[,1],poisson.density,lwd=2)
lines(foo[,1],nb.density,lwd=2,col="red")
lines(foo[,1],geo.density,lwd=2,col="blue")
legend("topright",lwd=2,col=c("black","red", "blue"),legend=c("Poisson","Negative Binomial", "Geometric"))

# Bins to zones
bindata <- data[, c('xbin', 'ybin')]
setDT(bindata)[, paste0("x", 1:2) := tstrsplit(xbin, split = ",", fixed = TRUE)]
setDT(bindata)[, paste0("y", 1:2) := tstrsplit(ybin, split = ",", fixed = TRUE)]
bindata[, x1 := as.numeric(substr(x1, 2,3))]
bindata[, x2 := as.numeric(substr(x2, 1,2))]
bindata[, y1 := as.numeric(substr(y1, 2,3))]
bindata[, y2 := as.numeric(substr(y2, 1,2))]
bindata[, xmid := 0.5*(x1 + x2)]
bindata[, ymid := 0.5*(y1 + y2)]
bindata <- bindata[, c('xbin', 'ybin', "xmid", "ymid") ]

bindata <- bindata[, .N, by = c('xbin', 'ybin', "xmid", "ymid")]
bindata <- bindata[bindata$N >= 8]
bindata$N <- NULL

setorder(bindata, -ybin, xbin)
bindata$zone <- 1:nrow(bindata)

data <- bindata[data, , on = c('xbin', 'ybin'), nomatch=NULL]
data <- data[, c('time','type', 'x', 'y', 'zone')]

save(data, type_meta_data, bindata, file = 'dataset.RData')
load("~/Documents/MyResearch/Terrorism/Code/dataset.RData")
source("~/Documents/MyResearch/Terrorism/Code/helpers.R")
#source("~/Documents/MyResearch/Terrorism/Code/helpers_bessell.R")

# # probability of atleast one attack per day
# events <- data$time[data$zone <= 13]
# events <- unique(events)
# N <- 2000
# p_hat <- length(events[events <= N])/N
# 
# # K function for atleast one event per day
# K_hats <- sapply(1:100, function(z) InHomoDiscK(s = z, events = events, N = N, p_hat = p_hat) - z)
# plot(K_hats)

# PHASE I: ONLY TIMES

# probability of atleast one attack per day
events <- data$time[data$zone <= 9]
N <- 2000
p_hat <- length(events[events <= N])/N

# K function
K_hats <- sapply(1:100, function(z) InHomoDiscKmulti(s = z, events = events, N = N, p_hat = p_hat) - z)
plot(K_hats)

# Simulated test on non-clustered data
k <- 2500
events <- rep(1:k, times = rpois(k, p_hat))
N <- 2000

K_hats <- sapply(1:100, function(z) InHomoDiscKmulti(s = z, events = events, N = N, p_hat = p_hat) - z)
plot(K_hats)

# Experiments

# Training data
N <- 2000
events = data$time[data$time <= N & data$zone <= 9]

# Baseline
poisson_rate = length(events)/N
vec <- rep(0, N)
names(vec) <- 1:N
vec[names(table(events))] <- unlist(table(events))
base_nb <- fitdistr(vec,"negative binomial")

# MLE
par = trans_par(c(0.1, 0.75, 0.05))
#par = trans_par(c(0.1, 0.75, 200))

estimates <-
  nlm(
    p = par,
    f = mle_loglik,
    events = events,
    fscale = mle_loglik(par, events),
    gradtol = 1e-6,
    steptol = 1e-6,
    print.level = 2,
    hessian = T
  )

mle_params = round(inv_trans_par(estimates$estimate), 4)

mle_std_errors = sqrt(diag(solve(estimates$hessian)))

lower = round(inv_trans_par(estimates$estimate - 1.96*mle_std_errors), 4)
upper = round(inv_trans_par(estimates$estimate + 1.96*mle_std_errors), 4)

rbind(lower, mle_params, upper)

# Params marginal gut check
events_sim = simulator(params = mle_params, N = max(data$time))
vec <- rep(0, max(events_sim))
names(vec) <- 1:max(events_sim)
vec[names(table(events_sim))] <- unlist(table(events_sim))
simulated_distribution = table(vec)/length(vec)

xx  <- 0:max(vec)
counts <- as.vector(table(factor(vec, levels = xx), exclude = NULL))
foo <- barplot(counts,names.arg=xx,ylim=range(counts))

# Validation
truth_vec <- vector()
lambda_vec <- vector()
for(i in (N+1):(N+922)){
  truth_i = length(data$time[data$time == i & data$zone <= 9])
  history_i = data$time[data$time < i & data$zone <= 9]
  lambda_vec = c(lambda_vec, simulator_history(params = mle_params, history = history_i, N = i-1))
  truth_vec = c(truth_vec, truth_i)
}

model_scores = dpois(x = truth_vec, lambda = lambda_vec, log = T)
poisson_scores = dpois(x = truth_vec, lambda = poisson_rate, log = T)
nb_scores = dnbinom(x = truth_vec, size=base_nb$estimate["size"], mu=base_nb$estimate["mu"], log = T)

sum(model_scores)
sum(poisson_scores)
sum(nb_scores)

# PHASE II: TIMES AND ZONES
load("~/Documents/MyResearch/Terrorism/Code/dataset.RData")
source("~/Documents/MyResearch/Terrorism/Code/helpers_stpp_simple.R")
library(mvtnorm)

# Experiments

# Training data
N = 2000
Z = 100
my_sigma = 0.01

extra_columns <- data.frame()
for(i in 1:nrow(data)){
  extra_columns <- rbind(extra_columns , sapply(1:Z, function(z) pmvnorm(lower = c(bindata$xmid[z] - 0.5, bindata$ymid[z] - 0.5), 
                                                                         upper = c(bindata$xmid[z] + 0.5, bindata$ymid[z] + 0.5), 
                                                                         mean = c(data$x[i], data$y[i]), sigma = my_sigma*diag(2))) )
}
colnames(extra_columns) <- 1:Z
extra_columns[extra_columns < 1e-6] <- 0
extra_columns = extra_columns/rowSums(extra_columns)
data <- cbind(data, extra_columns)

events = data[data$time <= N]

# Baseline
poisson_rate = as.numeric(table(events$zone))/N
nb_mat <- as.matrix.data.frame(table(events$time, events$zone))
nb_mat <- rbind(nb_mat, matrix(0, nrow = N - nrow(nb_mat), ncol = nrow(bindata)))
base_nb <- append(append(apply(nb_mat[, 1:12], 2, function(z) fitdistr(z,"negative binomial", start = list(size = 0.1, mu = 0.1))),
                  apply(nb_mat[, 13, drop = F], 2, function(z) fitdistr(z,"negative binomial", start = list(size = 30, mu = 0.9)))),
                apply(nb_mat[, -c(1:13)], 2, function(z) fitdistr(z,"negative binomial", start = list(size = 0.1, mu = 0.1))))

# MLE
par = trans_par( c( 0.0049, 0.8805, 0.0213) )

# Initialization
event_mat <- as.matrix(events[, lapply(.SD, sum), .SDcols = 6:(Z+5), by = 'time'])
event_mat <- rbind(event_mat, cbind(time = (1:N)[!1:N %in% event_mat[,1]], matrix(0, nrow = N - nrow(event_mat), ncol = Z)))
event_mat <- event_mat[order(event_mat[,1]),]
count_mat = as.data.frame.matrix(table(factor(events$time, levels = 1:N), events$zone))
#count_vec = rowSums(count_mat)
mle_loglik(par, event_mat, count_mat, N, Z)

estimates <-
  nlm(
    p = par,
    f = mle_loglik,
    event_mat = event_mat,
    count_mat = count_mat,
    N = N,
    Z = Z,
    fscale = mle_loglik(par, event_mat, count_mat, N, Z),
    gradtol = 1e-6,
    steptol = 1e-6,
    print.level = 2,
    hessian = T
  )

mle_params = round(inv_trans_par(estimates$estimate), 4)

mle_std_errors = sqrt(diag(solve(estimates$hessian)))

lower = round(inv_trans_par(estimates$estimate - 1.96*mle_std_errors), 4)
upper = round(inv_trans_par(estimates$estimate + 1.96*mle_std_errors), 4)

rbind(lower, mle_params, upper)

# Validation
data$zone = factor(data$zone, levels = 1:Z)
truth_mat  <- data.frame()
lambda_mat <- data.frame()

for(i in (N+1):(N+922)){
  truth_i = as.numeric(table(data$zone[data$time == i]))
  history_i = data[data$time < i]
  event_mat = as.matrix(history_i[, lapply(.SD, sum), .SDcols = 6:(Z+5), by = 'time'])
  event_mat = rbind(event_mat, cbind(time = (1:(i-1))[!1:(i-1) %in% event_mat[,1]], matrix(0, nrow = (i-1) - nrow(event_mat), ncol = Z)))
  event_mat = event_mat[order(event_mat[,1]),]
  count_mat  = as.data.frame.matrix(table(factor(history_i$time, levels = 1:(i-1)), history_i$zone))
  lambda_mat = rbind(lambda_mat, simulator_history(params = mle_params, event_mat, count_mat, N = i-1, Z))
  truth_mat  = rbind(truth_mat, truth_i)
}

model_scores = sapply(1:Z, FUN = function(z) dpois(truth_mat[, z], lambda_mat[, z], log = T))
poisson_scores = sapply(1:Z, FUN = function(z) dpois(truth_mat[, z], poisson_rate[z], log = T))
nb_scores = sapply(1:Z, FUN = function(z) dnbinom(truth_mat[, z], size=base_nb[[z]]$estimate["size"], 
                                                  mu=base_nb[[z]]$estimate["mu"],  log = T))

sum(model_scores)
sum(poisson_scores)
sum(nb_scores)

# PHASE III : TIMES, ZONES AND MARKS
load("~/Documents/MyResearch/Terrorism/Code/dataset.RData")
source("~/Documents/MyResearch/Terrorism/Code/helpers_mstpp.R")
library(mvtnorm)
library(doMC)
library(tscount)

# Experiments

# Training data
M = 6
N = 2000
Z = 100
my_sigma = 0.01

extra_columns <- data.frame()
for(i in 1:nrow(data)){
  extra_columns <- rbind(extra_columns , sapply(1:Z, function(z) pmvnorm(lower = c(bindata$xmid[z] - 0.5, bindata$ymid[z] - 0.5), 
                                                                         upper = c(bindata$xmid[z] + 0.5, bindata$ymid[z] + 0.5), 
                                                                         mean = c(data$x[i], data$y[i]), sigma = my_sigma*diag(2))) )
}
colnames(extra_columns) <- 1:Z
extra_columns[extra_columns < 1e-6] <- 0
extra_columns = extra_columns/rowSums(extra_columns)
data <- cbind(data, extra_columns)

events = data[data$time <= N]

# Baseline
count_array = unclass(table(factor(events$time, levels = 1:N), events$zone))
tsfit_pois_list <- rep(list(0), 100)
for(i in 1:Z){
   if(sum(count_array[,i]) < 50){
     tsfit_pois_list[[i]] <- tsglm(ts(count_array[,i]), model = list(past_obs = 1, past_mean = 1),
                                   distr = "poisson")
   }else{
     tsfit_pois_list[[i]] <- tsglm(ts(count_array[,i]), model = list(past_obs = 1:2, past_mean = 1:2),
                                   distr = "poisson")
   }
}

# MLE
par = trans_par( c( 0.0008, 0.8805, 0.0213, rep(0.166, 30)), M = M )

par <- c(-7.1521154,  2.4200398, -3.9131329,  
         3.0084328,  0.2282321,  2.2169403,  0.3996411,  0.3758740,
         0.9609641,  0.6214124,  0.9409980, -0.8799363, -7.2413019,
         2.5171575,  0.7375440,  4.6011011,  0.6009406,  2.1562250,
         0.7532594, -0.8296697,  1.2226052,  2.0212997, -0.1845317,
         1.7508704,  0.1626820,  2.1401129,  0.7241631,  3.1523136,
         1.5529734, -6.8732359,  1.6618246,  0.3684000,  1.0026676)

# Initialization
event_mat_list <- rep(list(0), 6)
for(m in 1:M){
  event_mat <- as.matrix(events[events$type == m, lapply(.SD, sum), .SDcols = 6:(Z+5), by = 'time'])
  event_mat <- rbind(event_mat, cbind(time = (1:N)[!1:N %in% event_mat[,1]], matrix(0, nrow = N - nrow(event_mat), ncol = Z)))
  event_mat <- event_mat[order(event_mat[,1]),]
  
  event_mat_list[[m]] <- event_mat
}
count_array = unclass(table(factor(events$time, levels = 1:N), events$zone, events$type))
registerDoMC(cores = 6)

system.time(mle_loglik(par, event_mat_list, count_array, N, Z, M))
system.time(mle_loglik_pll(par, event_mat_list, count_array, N, Z, M))

estimates <-
  nlm(
    p = par,
    f = mle_loglik_pll,
    event_mat_list = event_mat_list,
    count_array = count_array,
    N = N,
    Z = Z,
    M = M,
    fscale = mle_loglik_pll(par, event_mat_list, count_array, N, Z, M),
    gradtol = 1e-6,
    steptol = 1e-6,
    iterlim = 100,
    print.level = 2,
    hessian = T
  )

mle_params = round(inv_trans_par(estimates$estimate, M = M), 6)

mle_std_errors = sqrt(diag(solve(estimates$hessian)))

lower = round(inv_trans_par(estimates$estimate - 1.96*mle_std_errors, M = M), 6)
upper = round(inv_trans_par(estimates$estimate + 1.96*mle_std_errors, M = M), 6)

rbind(lower, mle_params, upper)

# Parameters

mle_params[1]
mle_params[2]
sum(1 - (dgeom(x = 1:100, prob = mle_params[3])*mle_params[2]*mean(apply(extra_columns, 1, max)) < mle_params[1]*M))

O_mark_mle  = matrix(mle_params[4:(M*(M-1)+3)], nrow = M, byrow = T)
O_mark_mle  = cbind(O_mark_mle, 1 - rowSums(O_mark_mle)) 
O_mark_mle

type_meta_data <- type_meta_data[order(type_meta_data$type),]
type_meta_data$short_type_label <- c('ARMED ASSAULT', 'ASSASSINATION', 'BOMBING', 'FACILITY ATTACK', 'KIDNAPPING', 'OTHERS')
row.names(O_mark_mle) <- colnames(O_mark_mle) <- type_meta_data$short_type_label
library(corrplot)
corrplot(O_mark_mle,
         method="circle",
         is.corr=FALSE,
         type="full",
         tl.srt = 45,
         number.cex = 1,
         cl.lim=c(0,1), 
         addCoef.col = rgb(0,0,0, alpha = 0.6)
)

# Validation
data$zone = factor(data$zone, levels = 1:Z)
truth_array  <- array(0, dim = c(922, Z, M))
lambda_array <- array(0, dim = c(922, Z, M))

event_mat_list <- rep(list(0), 6)

for(i in (N+1):(N+922)){
  truth_array[i-N,,] = unclass(table(data$zone[data$time == i], data$type[data$time == i]))
  history_i = data[data$time < i]
  for(m in 1:M){
    event_mat <- as.matrix(history_i[history_i$type == m, lapply(.SD, sum), .SDcols = 6:(Z+5), by = 'time'])
    event_mat <- rbind(event_mat, cbind(time = (1:(i-1))[!1:(i-1) %in% event_mat[,1]], matrix(0, nrow = (i-1) - nrow(event_mat), ncol = Z)))
    event_mat <- event_mat[order(event_mat[,1]),]
    
    event_mat_list[[m]] <- event_mat
  }
  count_array = unclass(table(factor(history_i$time, levels = 1:(i-1)), history_i$zone, history_i$type))
  lambda_array[i-N,,] = simulator_history(params = mle_params, event_mat_list, count_array, N = i-1, Z = Z, M = M)

}

# Baseline
baseline_array <- array(0, dim = c(922, Z))
count_array = unclass(table(factor(history_i$time, levels = 1:(i-1)), history_i$zone))
for(j in 1:Z){
    newobs = ts(count_array[,j])
    baseline_array[,j] = predict(tsfit_pois_list[[j]], n.ahead = i-N, newobs = newobs[(N+1):(i-1)], level = 0, global = TRUE)$pred
}

model_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], lambda_array[, z, m], log = T))))
mark_prop = unclass(table(events$type))/nrow(events)
baseline_scores = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], baseline_array[, z]*mark_prop[m], log = T))))
poisson_rate_3d = unclass(table(events$zone, events$type))/N
poisson_rate_3d[poisson_rate_3d < 1e-16] <- 1e-16
poisson_scores_3d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_3d[z, m], log = T))))
poisson_rate_2d = as.numeric(table(events$zone))/(N*M)
poisson_scores_2d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d[z], log = T))))
poisson_rate_1d = nrow(events)/(N*Z*M)
poisson_scores_1d = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_1d, log = T))))
poisson_rate_2d_clever = as.numeric(table(events$zone))/N
poisson_scores_2d_clever = sapply(1:M, FUN = function(m) sapply(1:Z, FUN = function(z) sum(dpois(truth_array[, z, m], poisson_rate_2d_clever[z]*mark_prop[m], log = T))))

sum(model_scores)
sum(baseline_scores)
sum(poisson_scores_2d_clever)
sum(poisson_scores_2d)
sum(poisson_scores_3d)
sum(poisson_scores_1d)

# Predictive model assessment

model_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                          function(z) scoring(response = truth_array[, z, m], pred = lambda_array[, z, m], 
                                                distr = "poisson"))))

baseline_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                                function(z) scoring(response = truth_array[, z, m], pred = baseline_array[, z]*mark_prop[m], 
                                                      distr = "poisson"))))

poisson_scores_array  = sapply(1:M, FUN = function(m) rowSums(sapply(1:Z, FUN = 
                                                 function(z) scoring(response = truth_array[, z, m], pred = rep(poisson_rate_2d_clever[z]*mark_prop[m], 922), 
                                                       distr = "poisson"))))

pred_perf <- as.data.frame(cbind(rowSums(model_scores_array), rowSums(baseline_scores_array),rowSums(poisson_scores_array)))
colnames(pred_perf) <- c('DTMSTPP', 'TSGLM', 'POISSON')
pred_perf$scoring <- row.names(pred_perf)

write.table(pred_perf[,c(4, 1:3)], file = "pred_perf.csv", row.names = F)

# Prediction Accuracy Index

# Hotspot mapping
dummy1 <- dgeom(1:100, prob = mle_params[3])
dummy2 <- 1/(50 + ((1:100)*3) )
dummy2/dummy1
cbind(dummy1, dummy2)

hotspot_array <- array(0, dim = c(922, Z))
N <- 2000
for(i in (N+1):(N+922)){
  history_i = data[data$time < i]
  time_weights = (1/(1 + (i - history_i$time)/7))
  time_weights[time_weights < 0.065] = 0
  for(j in 1:Z){
   dist_weights = 1/(1 + (abs(bindata$xmid[j] - history_i$x) %/% 0.5) +  (abs(bindata$ymid[j] - history_i$y) %/% 0.5) )
   dist_weights[dist_weights <= .2] = 0
   hotspot_array[i-N, j] = sum(time_weights*dist_weights)
  }
}

hotspot_pais <- array(0, dim = c(70, 2))
theta <- (seq(0, 0.7, 0.01)^3)*(round(max(hotspot_array)+1, 0))
for(i in 1:70){
    hotspot_pais[i,] <- hotspot_PAI(truth_array, hotspot_array, mark_prop, theta[i])
}
hotspot_pais <- data.frame(cbind(hotspot_pais, hotspot_pais[,1]/hotspot_pais[,2]))
colnames(hotspot_pais) <- c('hit_rate', 'area_perc', 'PAI')
plot(x = hotspot_pais$hit_rate, y = hotspot_pais$PAI)

# Validation
model_pais <- baseline_pais <- array(0, dim = c(100, 100, 2))
theta <- (seq(0, 1, 0.01)^5)*0.5
zeta <- (seq(0, 1, 0.01)^5)*0.05
for(i in 1:100){
  for(j in 1:100){
    model_pais[i,j,] <- model_PAI(truth_array, lambda_array, theta[i], zeta[j])
    baseline_pais[i,j,] <- baseline_PAI(truth_array, baseline_array, mark_prop, theta[i], zeta[j])
  }
}
m_pais_mat <- matrix(model_pais, 10000, 2)
b_pais_mat <- matrix(baseline_pais, 10000, 2)
plot(m_pais_mat[,2], m_pais_mat[,1])
plot(b_pais_mat[,2], b_pais_mat[,1])
m_pais_mat <- cbind(m_pais_mat, m_pais_mat[,1]/m_pais_mat[,2])
b_pais_mat <- cbind(b_pais_mat, b_pais_mat[,1]/b_pais_mat[,2])
plot(m_pais_mat[,1], m_pais_mat[,3])
plot(b_pais_mat[,1], b_pais_mat[,3])

# Unique vals
m_pais_mat <- data.table(m_pais_mat)
b_pais_mat <- data.table(b_pais_mat)
colnames(m_pais_mat) <- colnames(b_pais_mat) <- c('hit_rate', 'area_perc', 'PAI')
m_pais_mat_upd <- m_pais_mat[, .SD[which.max(PAI)], by = c('hit_rate','area_perc')]
b_pais_mat_upd <- b_pais_mat[, .SD[which.max(PAI)], by = c('hit_rate','area_perc')]
m_pais_mat_upd$hit_rate_bin <- cut(m_pais_mat_upd$hit_rate, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
m_pais_mat_upd$area_perc_bin <- cut(m_pais_mat_upd$area_perc, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
b_pais_mat_upd$hit_rate_bin <- cut(b_pais_mat_upd$hit_rate, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
b_pais_mat_upd$area_perc_bin <- cut(b_pais_mat_upd$area_perc, breaks = seq(0, 1, 0.02), labels = seq(0.01, 0.99, 0.02), include.lowest = T)
m_pais_mat_upd <- m_pais_mat_upd[, .SD[which.max(PAI)], by = c('hit_rate_bin')]
b_pais_mat_upd <- b_pais_mat_upd[, .SD[which.max(PAI)], by = c('hit_rate_bin')]

plot(x = m_pais_mat_upd$area_perc, y = m_pais_mat_upd$hit_rate)
plot(x = m_pais_mat_upd$area_perc, y = m_pais_mat_upd$PAI)
plot(x = m_pais_mat_upd$hit_rate, y = m_pais_mat_upd$PAI)

plot(x = b_pais_mat_upd$area_perc, y = b_pais_mat_upd$hit_rate)
plot(x = b_pais_mat_upd$area_perc, y = b_pais_mat_upd$PAI)
plot(x = b_pais_mat_upd$hit_rate,  y = b_pais_mat_upd$PAI)

# plot data
library(hrbrthemes)
plot_data <- rbind(cbind(m_pais_mat_upd[, c('hit_rate', 'PAI')], 'DTMSTPP'), cbind(b_pais_mat_upd[, c('hit_rate', 'PAI')], 'Hotspot Mapping'))
colnames(plot_data) <- c('HIT.RATE', 'PAI', 'MODEL')

ggplot(plot_data, aes(x=HIT.RATE, y=PAI, colour = MODEL)) +
  geom_line() + 
  xlab("Hit Rate") +
  ylab("Prediction Accuray Index") +
  theme_ipsum(base_size = 16,
              axis_title_size = 16) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'bottom')


# branching structure
branch_probs <- branching_structure(params = mle_params, events = as.data.frame(data), M = M)
mode_struct = apply(branch_probs, 1, FUN = function(z) which.max(z) - 1)
branch_data <- cbind(data[, 1:5], b_str = mode_struct)

long_chain <- function(n, branch_data){
  if(branch_data$b_str[n] == 0){
    return(0)
  }else{
    return(1 + long_chain(branch_data$b_str[n], branch_data)) 
  }
}

result_vec <- rep(0, 5923)
for(i in 1:5923){
  result_vec[i] <- long_chain(i, branch_data)
}

# Hotspot mapping