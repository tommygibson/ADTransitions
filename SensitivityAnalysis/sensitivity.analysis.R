###### Sensitivity analysis after grid search

##### Initial values for the grid search "winner" were (-12, 0.02) for (log(k0), k1)
##### Now, we redo optimization as a sensitivity analysis
##### Lower: multiply incidence rates by 0.5, multiply r45 by 0.66 and set r5,10 = 0.26
##### Higher: multiply incidence by 1.5, multiply r45 by 1.65, and set r5,10 = 0.34

library(tidyverse)
library(ggplot2)
library(knitr)
library(readxl)
library(nloptr)
library(here)


# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Brookmeyer/Code")
source(here("BrookFuncs.R"))
source(here("AD_eval.f.g.R"))

ages <- 50:95
inc.ages <- 65:90

empirical.incidence <- 0.00117 * exp(0.126 * (inc.ages - 60)) * 100
incidence.low <- empirical.incidence * 0.5
incidence.high <- empirical.incidence * 1.5

sens.inits <- rep(c(-12, 0.02), 12)


## lower bounds and upper bounds for parameters
lb <- rep(c(-14, 0.001), 12)

ub <- rep(c(-4, 0.15), 12)

opt.low <- nloptr(x0 = sens.inits,
                  eval_f = eval_f_logs_weighted_low,
                  lb = lb, ub = ub,
                  eval_g_ineq = eval_g_ineq_weighted,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1e-3,
                             "maxeval"=40000),
                  r45 = r45.params,
                  prevs = avg_prev_u,
                  incidence = incidence.low,
                  w = 1)

opt.high <- nloptr(x0 = sens.inits,
                   eval_f = eval_f_logs_weighted_high,
                   lb = lb, ub = ub,
                   eval_g_ineq = eval_g_ineq_weighted,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=40000),
                   r45 = r45.params,
                   prevs = avg_prev_u,
                   incidence = incidence.high,
                   w = 1)

opt.sensitivity <- list(opt.low, opt.high)
names(opt.sensitivity) <- c("opt.low", "opt.high")
saveRDS(opt.sensitivity, "OptResults/opt.sensitivity.rds")

#### Make lifetime tables so we don't need to keep recalculating
tab.ages <- seq(60, 90, 5)
low.mats <- make_trans_matrix_low(opt.low$solution, r45 = r45.params)
high.mats <- make_trans_matrix_high(opt.high$solution, r45 = r45.params)

lifetime.table.f.low <- as.data.frame(matrix(nrow = length(tab.ages), ncol = 10))
lifetime.table.m.low <- as.data.frame(matrix(nrow = length(tab.ages), ncol = 10))
lifetime.table.f.high <- as.data.frame(matrix(nrow = length(tab.ages), ncol = 10))
lifetime.table.m.high <- as.data.frame(matrix(nrow = length(tab.ages), ncol = 10))

lifetime.table.f.low[,1] <- lifetime.table.m.low[,1] <- 
  lifetime.table.f.high[,1] <- lifetime.table.m.high[,1] <- tab.ages


for(i in 1:length(tab.ages)){
  for(j in 1:9){
    
    curr.f.low <- lifetime(age = tab.ages[i], g = "Female", state = j, k0 = low.mats[[1]], k1 = low.mats[[2]])
    curr.m.low <- lifetime(age = tab.ages[i], g = "Male", state = j, k0 = low.mats[[1]], k1 = low.mats[[2]])
    curr.f.high <- lifetime(age = tab.ages[i], g = "Female", state = j, k0 = high.mats[[1]], k1 = high.mats[[2]])
    curr.m.high <- lifetime(age = tab.ages[i], g = "Male", state = j, k0 = high.mats[[1]], k1 = high.mats[[2]])
    
    lifetime.table.f.low[i, (j + 1)] <- curr.f.low[1]
    lifetime.table.m.low[i, (j + 1)] <- curr.m.low[1]
    lifetime.table.f.high[i, (j + 1)] <- curr.f.high[1]
    lifetime.table.m.high[i, (j + 1)] <- curr.m.high[1]
    
  }
}

sensitivity.lifetime.tables <- list(lifetime.table.f.low, lifetime.table.m.low, lifetime.table.f.high, lifetime.table.m.high)
names(sensitivity.lifetime.tables) <- c("f.low", "m.low", "f.high", "m.high")
for(i in 1:length(sensitivity.lifetime.tables)){
  names(sensitivity.lifetime.tables[[i]]) <- c("Age", "Normal", "A", "A+T", "A+T+N",
                                               "A+T+N + MCI", "T", "T+N", "N", "A+N")
}

saveRDS(sensitivity.lifetime.tables, "SensitivityAnalysis/sensitivity.lifetime.tables.rds")



