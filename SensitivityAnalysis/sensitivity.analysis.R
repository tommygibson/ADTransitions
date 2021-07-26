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
source(here("MCI.prevalence.R"))

ages <- 50:95
inc.ages <- 65:90

incidence <- 0.00117 * exp(0.126 * (inc.ages - 60)) * 100
incidence.low <- incidence * 0.5
incidence.high <- incidence * 1.5

sens.inits <- rep(c(-12, 0.02), 12)


## lower bounds and upper bounds for parameters
lb <- rep(c(-14, 0.001), 12)

ub <- rep(c(-4, 0.15), 12)

lb.test <- rep(c(-15, 0.00001), 12)
ub.test <- rep(c(-0.00001, 0.20), 12)

opt.middle.test <- nloptr(x0 = sens.inits,
                          eval_f = eval_f_logs_weighted,
                          lb = lb.test, ub = ub.test,
                          eval_g_ineq = eval_g_ineq_weighted,
                          opts = list("algorithm"="NLOPT_LN_COBYLA",
                                      "xtol_rel"=1e-3,
                                      "maxeval"=40000),
                          r45 = r45.params,
                          prevs = avg_prev_u,
                          incidence = incidence,
                          w = 1)

opt.low <- nloptr(x0 = sens.inits,
                  eval_f = eval_f_logs_weighted_low,
                  lb = lb.test, ub = ub.test,
                  eval_g_ineq = eval_g_ineq_weighted,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1e-3,
                             "maxeval"=40000),
                  r45 = r45.params,
                  prevs = avg_prev_u_low,
                  incidence = incidence.low,
                  w = 1)

opt.high <- nloptr(x0 = sens.inits,
                   eval_f = eval_f_logs_weighted_high,
                   lb = lb.test, ub = ub.test,
                   eval_g_ineq = eval_g_ineq_weighted,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=40000),
                   r45 = r45.params,
                   prevs = avg_prev_u_high,
                   incidence = incidence.high,
                   w = 1)

opt.sensitivity <- list(opt.middle.test, opt.low, opt.high)
names(opt.sensitivity) <- c("opt.middle", "opt.low", "opt.high")
saveRDS(opt.sensitivity, "OptResults/opt.sensitivity.unrestricted.rds")

#### Make lifetime tables so we don't need to keep recalculating
tab.ages <- seq(60, 90, 5)
r45.low <- r45.high <- r45.params
r45.low[1] <- r45.params[1] * 0.66
r45.high[1] <- r45.params[1] * 1.65
low.mats <- make_trans_matrix_low(opt.low$solution, r45 = r45.low)
mid.mats <- make_trans_matrix(opt.middle.test$solution, r45 = r45.params)
high.mats <- make_trans_matrix_high(opt.high$solution, r45 = r45.high)

lifetime.table.f.low <- lifetime.table.m.low <- 
  lifetime.table.f.mid <- lifetime.table.m.mid <- 
  lifetime.table.f.high <- lifetime.table.m.high <- as.data.frame(matrix(nrow = length(tab.ages), ncol = 10))

lifetime.table.f.low[,1] <- lifetime.table.m.low[,1] <- 
  lifetime.table.f.high[,1] <- lifetime.table.m.high[,1] <- tab.ages


for(i in 1:length(tab.ages)){
  for(j in 1:9){
    
    curr.f.low <- lifetime(age = tab.ages[i], g = "Female", state = j, k0 = low.mats[[1]], k1 = low.mats[[2]])
    curr.m.low <- lifetime(age = tab.ages[i], g = "Male", state = j, k0 = low.mats[[1]], k1 = low.mats[[2]])
    curr.f.mid <- lifetime(age = tab.ages[i], g = "Female", state = j, k0 = mid.mats[[1]], k1 = mid.mats[[2]])
    curr.m.mid <- lifetime(age = tab.ages[i], g = "Male", state = j, k0 = mid.mats[[1]], k1 = mid.mats[[2]])
    curr.f.high <- lifetime(age = tab.ages[i], g = "Female", state = j, k0 = high.mats[[1]], k1 = high.mats[[2]])
    curr.m.high <- lifetime(age = tab.ages[i], g = "Male", state = j, k0 = high.mats[[1]], k1 = high.mats[[2]])
    
    lifetime.table.f.low[i, (j + 1)] <- curr.f.low[1]
    lifetime.table.m.low[i, (j + 1)] <- curr.m.low[1]
    lifetime.table.f.mid[i, (j + 1)] <- curr.f.mid[1]
    lifetime.table.m.mid[i, (j + 1)] <- curr.m.mid[1]
    lifetime.table.f.high[i, (j + 1)] <- curr.f.high[1]
    lifetime.table.m.high[i, (j + 1)] <- curr.m.high[1]
    
  }
}

sensitivity.lifetime.tables <- list(lifetime.table.f.low, lifetime.table.m.low, 
                                    lifetime.table.f.mid, lifetime.table.m.mid,
                                    lifetime.table.f.high, lifetime.table.m.high)

names(sensitivity.lifetime.tables) <- c("f.low", "m.low", "f.mid", "m.mid", "f.high", "m.high")
for(i in 1:length(sensitivity.lifetime.tables)){
  names(sensitivity.lifetime.tables[[i]]) <- c("Age", "Normal", "A", "A+T", "A+T+N",
                                               "A+T+N + MCI", "T", "T+N", "N", "A+N")
}

saveRDS(sensitivity.lifetime.tables, "SensitivityAnalysis/sensitivity.lifetime.tables.rds")



