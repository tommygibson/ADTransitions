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

## Grid of initial values
Lk0.init <- seq(-12, -6, 2)
a <- length(Lk0.init)
k1.init <- seq(0.02, 0.1, length.out = a)

possible.inits <- as.matrix(expand.grid(Lk0.init, k1.init))


init.valid <- vector(length = (a * a))
init.loss <- vector(length = (a * a))

for(i in 1:(a * a)){
  g <- eval_g_ineq_both(rep(possible.inits[i,], 12), 0, 0, 0) < 0
  if(any(g == FALSE)){
    init.valid[i] <- 0
    init.loss[i] <- NA
  } else {
    init.valid[i] <- 1
    init.loss[i] <- eval_f_logs_weighted(x = rep(possible.inits[i,], 12), 
                                         r45.params, avg_prev_u, incidence = empirical.incidence, w = 1)
  }
}

###### the cobyla algorithm will sometimes take the parameter values outside
###### of the constraints, so this little precursor optimization accounts for that
for(i in 1:(a * a)){
  
  if(init.valid[i] == 0) next
  else{
    opt.temp <- nloptr(x0 = rep(possible.inits[i,], 12), 
                        eval_f = eval_f_logs_weighted, 
                        lb = lb, ub = ub, 
                        eval_g_ineq = eval_g_ineq_weighted,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1e-3,
                                    "maxeval"=10),
                        r45 = r45.params,
                        prevs = avg_prev_u,
                        incidence = empirical.incidence,
                        w = 1)
    if(is.nan(opt.temp$objective)){ 
      init.valid[i] <- 0
      init.loss[i] <- NA
    }
  }
  
}

opts.grid <- list()

index <- 1
for(i in 1:(a * a)){
  if(init.valid[i] == 0) next
  else{
    opts.grid[[index]] <- nloptr(x0 = rep(possible.inits[i,], 12), 
                               eval_f = eval_f_logs_weighted, 
                               lb = lb, ub = ub, 
                               eval_g_ineq = eval_g_ineq_weighted,
                               opts = list("algorithm"="NLOPT_LN_COBYLA",
                                           "xtol_rel"=5e-3,
                                           "maxeval"=15000),
                               r45 = r45.params,
                               prevs = avg_prev_u,
                               incidence = empirical.incidence,
                               w = 1)
    index <- index + 1
  }
}

#### Okay, we'll start with optimizing for 500 iterations, weed some out, 
# then 1000 iterations, weed out again, yada yada
##### List of all optimizations
# opts.500 <- list()
# opts.1000 <- list()
# index <- 1
# for(i in 1:(a * a)){
#   
#   if(init.valid[i] == 0) next
#   else{
#     opts.500[[index]] <- nloptr(x0 = rep(possible.inits[i,], 12), 
#                        eval_f = eval_f_logs_weighted, 
#                        lb = lb, ub = ub, 
#                        eval_g_ineq = eval_g_ineq_weighted,
#                        opts = list("algorithm"="NLOPT_LN_COBYLA",
#                                    "xtol_rel"=5e-3,
#                                    "maxeval"=500),
#                        r45 = r45.params,
#                        prevs = avg_prev_u,
#                        incidence = empirical.incidence,
#                        w = 1)
#     opts.1000[[index]] <- nloptr(x0 = rep(possible.inits[i,], 12), 
#                                 eval_f = eval_f_logs_weighted, 
#                                 lb = lb, ub = ub, 
#                                 eval_g_ineq = eval_g_ineq_weighted,
#                                 opts = list("algorithm"="NLOPT_LN_COBYLA",
#                                             "xtol_rel"=5e-3,
#                                             "maxeval"=1000),
#                                 r45 = r45.params,
#                                 prevs = avg_prev_u,
#                                 incidence = empirical.incidence,
#                                 w = 1)
#     index <- index + 1
#   }
#   
# }
# 
# ##### Do it again with 3000 iterations, only those qualifying are in
# opts.2000 <- list()
# 
# index <- 1
# for(i in 1:(a * a)){
#   if(init.valid[i] == 0) next
#   else{
#     opts.2000[[index]] <- nloptr(x0 = rep(possible.inits[i,], 12), 
#                                 eval_f = eval_f_logs_weighted, 
#                                 lb = lb, ub = ub, 
#                                 eval_g_ineq = eval_g_ineq_weighted,
#                                 opts = list("algorithm"="NLOPT_LN_COBYLA",
#                                             "xtol_rel"=1e-3,
#                                             "maxeval"=2000),
#                                 r45 = r45.params,
#                                 prevs = avg_prev_u,
#                                 incidence = empirical.incidence,
#                                 w = 1)
#     index <- index + 1
#   }
# }

####### save results
# saveRDS(opts.500, "OptResults/grid.500.rds")
# saveRDS(opts.1000, "OptResults/grid.1000.rds")
# saveRDS(opts.2000, "OptResults/grid.2000.rds")
#######

##### OBJECTIVE FUNCTION VALUES

obj.500 <- obj.1000 <- obj.2000 <- vector(length = length(opts.500))


for(i in 1:length(opts.500)){
  obj.500[i] <- opts.500[[i]]$objective
  obj.1000[i] <- opts.1000[[i]]$objective
  obj.2000[i] <- opts.2000[[i]]$objective
}

## make a pretty plot
obj.prog <- cbind.data.frame(c(init.loss[!is.na(init.loss)], obj.500, obj.1000, obj.2000),
                             rep(c(0, 500, 1000, 2000), each = sum(!is.na(init.loss))),
                             rep(1:sum(init.valid), 4))

names(obj.prog) <- c("Objective", "Iterations", "Init Set")

obj.prog %>%
  ggplot(aes(x = Iterations, y = Objective, group = `Init Set`)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_y_continuous(trans = "log2") +
  labs(title = "Convergence of COBYLA algorithm for multiple sets of initial values",
       y = "Objective function (log scale)") +
  scale_x_continuous(breaks = c(0, 500, 1000, 2000))


