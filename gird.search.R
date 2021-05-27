######## Grid search for global optimum


library(tidyverse)
library(ggplot2)
library(knitr)
library(readxl)
library(gridExtra)
library(grid)
library(gghighlight)
library(nloptr)


# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Brookmeyer/Code")
source("BrookFuncs.R")
source("AD_eval.f.g.R")

ages <- 50:95
inc.ages <- 65:90

empirical.incidence <- 0.00117 * exp(0.126 * (inc.ages - 60)) * 100


## lower bounds and upper bounds for parameters
lb <- rep(c(-14, 0.001), 12)

ub <- rep(c(-4, 0.15), 12)

Lk0.init <- seq(-12, -6, 2)
a <- length(Lk0.init)
k1.init <- seq(0.02, 0.1, length.out = a)

possible.inits <- as.matrix(expand.grid(Lk0.init, k1.init))

init.valid <- vector(length = (a * a))
init.loss <- vector(0, length = (a * a))

for(i in 1:(a * a)){
  g <- eval_g_ineq_both(rep(possible.inits[i,], 12), 0, 0, 0) < 0
  if(any(g == FALSE)){
    init.valid[i] <- 0
  } else {
    init.valid[i] <- 1
    init.loss[i] <- eval_f_logs_weighted(x = rep(possible.inits[i,], 12), 
                                      r45.params, avg_prev_u, incidence = empirical.incidence, w = 1)
  }
}

###### the cobyla algorithm will sometimes take the parameter values outside
###### of the constraints, so this little precursor optimization accounts for that
for(i in 1:(a * a)){
  
  if(init.valid[i] == 0 ) next
  else{
    opt.temp <- nloptr(x0 = rep(possible.inits[i,], 12), 
                        eval_f = eval_f_logs_weighted, 
                        lb = lb, ub = ub, 
                        eval_g_ineq = eval_g_ineq_weighted,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1e-3,
                                    "maxeval"=5),
                        r45 = r45.params,
                        prevs = avg_prev_u,
                        incidence = empirical.incidence,
                        w = 1)
    if(is.nan(opt.temp$objective)) init.valid[i] <- 0
    
  }
  
}

opts <- list()
index <- 1
for(i in 1:(a * a)){
  
  if(init.valid[i] == 0) next
  else{
    opts[[index]] <- nloptr(x0 = rep(possible.inits[i,], 12), 
                       eval_f = eval_f_logs_weighted, 
                       lb = lb, ub = ub, 
                       eval_g_ineq = eval_g_ineq_weighted,
                       opts = list("algorithm"="NLOPT_LN_COBYLA",
                                   "xtol_rel"=1e-3,
                                   "maxeval"=1000),
                       r45 = r45.params,
                       prevs = avg_prev_u,
                       incidence = empirical.incidence,
                       w = 1)
    index <- index + 1
  }
  
}

init.value.valid <- cbind(possible.inits, init.valid)

init.value.valid[,1:2][init.value.valid[,3] == 1]



