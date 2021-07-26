##### Optimization with reduced model

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
lb <- rep(c(-15, 0.0001), 9)

ub <- rep(c(-0.0001, 0.20), 9)

avg_prev_u_reduced <- as.matrix(
  as.data.frame(avg_prev_u) %>%
  mutate(women.nonA = `women.A-T+N-` + `women.A-T-N+` + `women.A-T+N+`) %>%
  select(-c("women.A-T+N-", "women.A-T-N+", "women.A-T+N+")) %>%
  relocate(women.nonA, .before = `women.A+T-N+`)
  )


opt.reduced <- nloptr(x0 = rep(c(-12, 0.02), 9), 
                      eval_f = eval_f_logs_weighted_reduced, 
                      lb = lb, ub = ub, 
                      eval_g_ineq = eval_g_ineq_weighted_reduced,
                      opts = list("algorithm"="NLOPT_LN_COBYLA",
                                  "xtol_rel"=1e-3,
                                  "maxeval"=30000),
                      r45 = r45.params,
                      prevs = avg_prev_u_reduced,
                      incidence = empirical.incidence,
                      w = 1)

saveRDS(opt.reduced, 'OptResults/opt.reduced.rds')

eval_f_logs_weighted_reduced(opt.reduced$solution, r45.params, avg_prev_u_reduced, empirical.incidence, w = 1)

k.red <- make_trans_matrix_reduced(opt.reduced$solution, r45.params)

prevs.red <- Prevrate.f.multi.ATN.uncond.reduced(age = 50:95, k0 = k.red[[1]], k1 = k.red[[2]])


prev.red <- cbind.data.frame(c(as.vector(prevs.red), as.vector(avg_prev_u_reduced)))


prev.red <- make_prevplot_data_f_reduced(opt.reduced$solution, r45.params = r45.params, prevs = avg_prev_u_reduced)
inc.red <- make_incplot_data_f_reduced(opt.reduced$solution, r45.params, empirical.incidence)

inc.red %>%
  ggplot(aes(x = Age, y = Incidence, color = Source)) +
  geom_line() +
  labs(title = "Empirical vs reduced multistate incidence rates",
       y = "Incidence (% per year)")

prev.red %>%
  ggplot(aes(x = Age, y = Prevalence, color = Source)) +
  geom_line() +
  facet_wrap(~State) +
  labs(title = "Reduced multistate prevalence rates",
       y = "Unconditional prevalence (proportion)")

