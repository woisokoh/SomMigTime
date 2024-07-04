burstiness_GB <- function(x){
  itv <- which(x > 0,arr.ind=TRUE)[,1]
  itv_dist <- diff(itv)
  r <- sd(itv_dist) / mean(itv_dist)
  B <- (r-1) / (r+1)
  
  n_tau <- length(itv_dist)
  tau_i1 <- itv_dist[1:(n_tau-1)]
  tau_i2 <- itv_dist[2:n_tau]
  m1 <- mean(tau_i1)
  m2 <- mean(tau_i2)
  sig1 <- sd(tau_i1)
  sig2 <- sd(tau_i2)
  M <- 0
  for (i in 1:(n_tau-1)){M = M +(tau_i1[i]-m1)*(tau_i2[i]-m2)}
  M <- M / (n_tau - 1)/ sig1 / sig2
  return(c(B,M))
}

pl_exp <- function(x){
  itv <- which(x > 0,arr.ind=TRUE)[,1]
  itv_dist <- diff(itv)
  
  pl <- displ$new(itv_dist)
  #est <- estimate_xmin(pl)
  pl$setXmin(1)
  est <- estimate_pars(pl)
  pl$setPars(est)
  
  expo <- disexp$new(itv_dist)
  #est <- estimate_xmin(expo)
  expo$setXmin(1)
  est <- estimate_pars(expo)
  expo$setPars(est)
  
  comp <- compare_distributions(pl,expo)
  return(c(comp$test_statistic, comp$p_two_sided))
}