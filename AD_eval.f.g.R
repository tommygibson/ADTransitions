### Holds functions for optimization in nloptr
### 
# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Brookmeyer/Code")
library(here)
source(here("BrookFuncs.R"))

eval_f_prev <- function(x, r45, prevs){
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  
  prevs.compare <- Prevrate.f.multi.ATN.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                               k0 = k0, k1 = k1)
  
  sumsquare <- sum((prevs.compare - prevs)^2)
  
  return(sumsquare)
  
}

# incidence loss function
eval_f_inc <- function(x, r45, incidence){
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  
  inc.compare <- incidence.f.multi(age = 65:90, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                   k0 = k0, k1 = k1) * 100
  
  sumsquare <- sum((log(inc.compare) - log(incidence))^2)
  
  return(sumsquare)
  
}

# both prevalence and incidence loss function
# no fancy weighting or anything, just straight up
# dont even use the same ages (50-95 prev, 65-90 incidence)

eval_f_both <- function(x, r45, prevs, incidence){
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- sum((prev.inc[, 1:8] - prevs) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

eval_f_both_logs <- function(x, r45, prevs, incidence){
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- 1 / (46 * 8) * sum((log(prev.inc[, 1:8] + .0001) - log(prevs + .0001)) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- (1 / 26) * sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

# w is the weight given to incidence, i.e. w = 1.5 means incidence is weighted 1.5 times as much as prev
eval_f_logs_weighted <- function(x, r45, prevs, incidence, w){
  
  w1 <- 1 / (1 + w)
  w2 <- w / (1 + w)
  
  n1 <- length(prevs)
  n2 <- length(incidence) 
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- w1 * (1 / n1) * sum((log(prev.inc[, 1:8] + 1e-4) - log(prevs + 1e-4)) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- w2 * (1 / n2) * sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

eval_f_logs_weighted_high <- function(x, r45, prevs, incidence, w){
  
  w1 <- 1 / (1 + w)
  w2 <- w / (1 + w)
  
  n1 <- length(prevs)
  n2 <- length(incidence) 
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- 1.65 * r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.34
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- w1 * (1 / n1) * sum((log(prev.inc[, 1:8] + 1e-4) - log(prevs + 1e-4)) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- w2 * (1 / n2) * sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

eval_f_logs_weighted_low <- function(x, r45, prevs, incidence, w){
  
  w1 <- 1 / (1 + w)
  w2 <- w / (1 + w)
  
  n1 <- length(prevs)
  n2 <- length(incidence) 
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- 0.66 * r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.26
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:95, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- w1 * (1 / n1) * sum((log(prev.inc[, 1:8] + 1e-4) - log(prevs + 1e-4)) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- w2 * (1 / n2) * sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

eval_f_logs_weighted_90 <- function(x, r45, prevs, incidence, w){
  
  w1 <- 1 / (1 + w)
  w2 <- w / (1 + w)
  
  n1 <- length(prevs)
  n2 <- length(incidence) 
  
  k0 <- matrix(0, nrow = 9, ncol = 10)
  k1 <- matrix(0, nrow = 9, ncol = 10)
  
  
  x[seq(1, 23, 2)] <- exp(x[seq(1, 23, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  k0[1, 8] <- x[5]
  k1[1, 8] <- x[6]
  # state 2
  k0[2, 3] <- x[7]
  k1[2, 3] <- x[8]
  k0[2, 9] <- x[9]
  k1[2, 9] <- x[10]
  # state 3
  k0[3, 4] <- x[11]
  k1[3, 4] <- x[12]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 10] <- 0.3
  k1[5, 10] <- 0
  # state 6
  k0[6, 3] <- x[13]
  k1[6, 3] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  # state 8
  k0[8, 7] <- x[19]
  k1[8, 7] <- x[20]
  k0[8, 9] <- x[21]
  k1[8, 9] <- x[22]
  # state 9
  k0[9, 4] <- x[23]
  k1[9, 4] <- x[24]
  
  prev.inc <- incidence.prevrate.f.uncond(age = 50:90, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, 
                                          k0 = k0, k1 = k1)
  
  # sums of squares for prevalences
  sumsquare.prev <- w1 * (1 / n1) * sum((log(prev.inc[, 1:8] + 1e-4) - log(prevs + 1e-4)) ^ 2)
  # sums of squares for log(incidence), only ages 65:90
  sumsquare.inc <- w2 * (1 / n2) * sum((log(prev.inc[16:41, 9] * 100) - log(incidence)) ^ 2)
  
  return(sumsquare.prev + sumsquare.inc)
  
}

###### Inequality constraints (each element of "constraints" should be < 0)
# last four constraints are on sum of transitions out of 1, 2, 6, 8

eval_g_ineq_prev <- function(x, r45, prevs){
  constraints <- c(x[1] + x[2] * 95,
                   x[3] + x[4] * 95,
                   x[5] + x[6] * 95,
                   x[7] + x[8] * 95,
                   x[9] + x[10] * 95,
                   x[11] + x[12] * 95,
                   x[13] + x[14] * 95,
                   x[15] + x[16] * 95,
                   x[17] + x[18] * 95,
                   x[19] + x[20] * 95,
                   x[21] + x[22] * 95,
                   x[23] + x[24] * 95,
                   log(exp(x[1] + x[2] * 95) + exp(x[3] + x[4] * 95) + exp(x[5] + x[6] * 95)),
                   log(exp(x[7] + x[8] * 95) + exp(x[9] + x[10] * 95)),
                   log(exp(x[13] + x[14] * 95) + exp(x[15] + x[16] * 95)),
                   log(exp(x[19] + x[20] * 95) + exp(x[21] + x[22] * 95)))
  return(constraints)
}

eval_g_ineq_inc <- function(x, r45, incidence){
  constraints <- c(x[1] + x[2] * 95,
                   x[3] + x[4] * 95,
                   x[5] + x[6] * 95,
                   x[7] + x[8] * 95,
                   x[9] + x[10] * 95,
                   x[11] + x[12] * 95,
                   x[13] + x[14] * 95,
                   x[15] + x[16] * 95,
                   x[17] + x[18] * 95,
                   x[19] + x[20] * 95,
                   x[21] + x[22] * 95,
                   x[23] + x[24] * 95,
                   log(exp(x[1] + x[2] * 95) + exp(x[3] + x[4] * 95) + exp(x[5] + x[6] * 95)),
                   log(exp(x[7] + x[8] * 95) + exp(x[9] + x[10] * 95)),
                   log(exp(x[13] + x[14] * 95) + exp(x[15] + x[16] * 95)),
                   log(exp(x[19] + x[20] * 95) + exp(x[21] + x[22] * 95)))
  return(constraints)
}

eval_g_ineq_both <- function(x, r45, prevs, incidence){
  constraints <- c(x[1] + x[2] * 95,
                   x[3] + x[4] * 95,
                   x[5] + x[6] * 95,
                   x[7] + x[8] * 95,
                   x[9] + x[10] * 95,
                   x[11] + x[12] * 95,
                   x[13] + x[14] * 95,
                   x[15] + x[16] * 95,
                   x[17] + x[18] * 95,
                   x[19] + x[20] * 95,
                   x[21] + x[22] * 95,
                   x[23] + x[24] * 95,
                   log(exp(x[1] + x[2] * 95) + exp(x[3] + x[4] * 95) + exp(x[5] + x[6] * 95)),
                   log(exp(x[7] + x[8] * 95) + exp(x[9] + x[10] * 95)),
                   log(exp(x[13] + x[14] * 95) + exp(x[15] + x[16] * 95)),
                   log(exp(x[19] + x[20] * 95) + exp(x[21] + x[22] * 95)))
  return(constraints)
}

eval_g_ineq_weighted <- function(x, r45, prevs, incidence, w){
  constraints <- c(x[1] + x[2] * 95,
                   x[3] + x[4] * 95,
                   x[5] + x[6] * 95,
                   x[7] + x[8] * 95,
                   x[9] + x[10] * 95,
                   x[11] + x[12] * 95,
                   x[13] + x[14] * 95,
                   x[15] + x[16] * 95,
                   x[17] + x[18] * 95,
                   x[19] + x[20] * 95,
                   x[21] + x[22] * 95,
                   x[23] + x[24] * 95,
                   log(exp(x[1] + x[2] * 95) + exp(x[3] + x[4] * 95) + exp(x[5] + x[6] * 95)),
                   log(exp(x[7] + x[8] * 95) + exp(x[9] + x[10] * 95)),
                   log(exp(x[13] + x[14] * 95) + exp(x[15] + x[16] * 95)),
                   log(exp(x[19] + x[20] * 95) + exp(x[21] + x[22] * 95)))
  return(constraints)
}


### Functions for creating plots of prevalence and incidence vs data
prev.ages <- 50:95
inc.ages <- 65:90

r45.low <- r45.high <- r45.params
r45.low[1] <- r45.params[1] * 0.66
r45.high[1] <- r45.params[1] * 1.65

make_prevplot_data_f <- function(sol, r45.params = r45.params, prevs = fem_prev_u){
  mats <- make_trans_matrix(sol, r45.params)
  
  prevs.u <- Prevrate.f.multi.ATN.uncond(prev.ages, k0 = mats[[1]], k1 = mats[[2]])
  
  dat.prev <- cbind.data.frame(c(as.vector(prevs.u), as.vector(prevs)),
                               rep(prev.ages, 8 * 2),
                               rep(rep(c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                         "T+", "T+ N+", "N+", "A+ N+"), each = length(prev.ages)), 2),
                               rep(c("Multistate", "Jack"), each = length(prev.ages) * 8))
  
  names(dat.prev) <- c("Prevalence", "Age", "State", "Source")
  
  dat.prev$State <- factor(dat.prev$State, levels = c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                      "T+", "T+ N+", "N+", "A+ N+"))
  
  return(dat.prev)
}

make_incplot_data_f <- function(sol, r45.params = r45.params, incidence.target = empirical.incidence){
  mats <- make_trans_matrix(sol, r45.params)
  
  inc <- incidence.f.multi(inc.ages, k0 = mats[[1]], k1 = mats[[2]]) * 100
  
  dat.inc <- cbind.data.frame(c(inc, incidence.target),
                              rep(inc.ages, 2),
                              rep(c("Multistate", "Empirical"), each = length(inc.ages)))
  
  names(dat.inc) <- c("Incidence", "Age", "Source")
  
  return(dat.inc)
}

make_prevplot_data_m <- function(sol, r45.params = r45.params, prevs = men_prev_u){
  mats <- make_trans_matrix(sol, r45.params)

  prevs.u <- Prevrate.m.multi.ATN.uncond(prev.ages, k0 = mats[[1]], k1 = mats[[2]])
  
  dat.prev <- cbind.data.frame(c(as.vector(prevs.u), as.vector(prevs)),
                               rep(prev.ages, 8 * 2),
                               rep(rep(c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                         "T+", "T+ N+", "N+", "A+ N+"), each = length(prev.ages)), 2),
                               rep(c("Multistate", "Jack"), each = length(prev.ages) * 8))
  
  names(dat.prev) <- c("Prevalence", "Age", "State", "Source")
  
  dat.prev$State <- factor(dat.prev$State, levels = c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                      "T+", "T+ N+", "N+", "A+ N+"))
  
  return(dat.prev)
}

make_incplot_data_m <- function(sol, r45.params = r45.params, incidence.target = empirical.incidence){
  mats <- make_trans_matrix(sol, r45.params)
  
  inc <- incidence.m.multi(inc.ages, k0 = mats[[1]], k1 = mats[[2]]) * 100
  
  dat.inc <- cbind.data.frame(c(inc, incidence.target),
                              rep(inc.ages, 2),
                              rep(c("Multistate", "Empirical"), each = length(inc.ages)))
  
  names(dat.inc) <- c("Incidence", "Age", "Source")
  
  return(dat.inc)
}

make_transitions <- function(Lk0_vec, k1_vec){
  ages <- 50:95
  transitions <- matrix(nrow = length(ages), ncol = length(Lk0_vec))
  
  for(i in 1:length(Lk0_vec)){
    transitions[,i] <- exp(Lk0_vec[i] + ages * k1_vec[i])
  }
  
  return(transitions)
  
}


