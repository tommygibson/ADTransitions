######## Optimization of AD dementia transition rate parameters using the nloptr package

library(tidyverse)
library(ggplot2)
library(knitr)

library(magrittr)
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




####### lower and upper bounds on parameters
# -20 < log(k0) < -0.1
# .005 < k1 < 0.4
lb <- rep(c(-14, 0.001), 12)

ub <- rep(c(-4, 0.15), 12)

###### Test constraints

# first set of inits: from optimization using demographic eqns
# giving base inits of k0ij = .0008 and k1 = .06 satisfies all constraints

init.simple <- rep(c(log(.0008), .06), 12)

eval_g_ineq_both(rep(c(log(.0008), .06), 12), r45.params, fem_prev_u, empirical.incidence) < 0 ## SUCCESS


##### Baseline loss from initial values

# base.both.loss <- eval_f_both(init.simple, r45.params, fem_prev_u, empirical.incidence)
# base.both.loss.logs <- eval_f_both_logs(init.simple, r45.params, fem_prev_u, empirical.incidence)
# # how much from prev vs incidence?
# # NOTE: not perfect match because incidence estimates will be from years prior to specified 'y'
# base.prev.loss <- eval_f_prev(init.simple, r45.params, fem_prev_u)
# base.inc.loss <- eval_f_inc(init.simple, r45.params, empirical.incidence)

# optimize just prevalence loss
# opt.prev <- nloptr(x0 = inits, 
#                eval_f = eval_f, 
#                lb = lb, ub = ub, 
#                eval_g_ineq = eval_g_ineq,
#                opts = list("algorithm"="NLOPT_LN_COBYLA",
#                            "xtol_rel"=1e-3,
#                            "maxeval"=5000),
#                r45 = r45.params,
#                prevs = fem_prev_u)
# 
# # optimize just incidence loss
# opt.inc <- nloptr(x0 = inits, 
#                   eval_f = eval_f_inc, 
#                   lb = lb, ub = ub, 
#                   eval_g_ineq = eval_g_ineq_inc,
#                   opts = list("algorithm"="NLOPT_LN_COBYLA",
#                               "xtol_rel"=1e-3,
#                               "maxeval"=5000),
#                   r45 = r45.params,
#                   incidence = empirical.incidence)

# optimize both prev and incidence
opt.both <- nloptr(x0 = init.simple,
                   eval_f = eval_f_both,
                   lb = lb, ub = ub,
                   eval_g_ineq = eval_g_ineq_both,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=15000),
                   r45 = r45.params,
                   prevs = fem_prev_u,
                   incidence = empirical.incidence)






###### k0 and k1 from optimization with both prev and incidence

opt.both.matrices <- make_trans_matrix(opt.both$solution, r45.params)

# incidence, unconditional prevalence, and conditional prevalence estimates using optimized parameters
opt.both.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014,
                                        k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]]) * 100
opt.both.prevalence.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
                                                   k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])
opt.both.prevalence.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1,
                                              k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])

# datasets for comparing with empirical incidence, jack prevalence
opt.both.inctest <- cbind.data.frame(c(opt.both.incidence, empirical.incidence),
                                     rep(inc.ages, 2),
                                     c(rep("Dual-optimized", length(inc.ages)),
                                       rep("Empirical", length(inc.ages))))
opt.both.prevtest <- cbind.data.frame(c(as.vector(opt.both.prevalence.u), as.vector(fem_prev_u)),
                                      c(as.vector(opt.both.prevalence.c), as.vector(fem_prev)),
                                      rep(ages, 16),
                                      rep(rep(c("State 1", "State 2", "State 3", "State 4",
                                                "State 6", "State 7", "State 8", "State 9"), each = length(ages)), 2),
                                      rep(c("Dual-optimized", "Jack"), each = length(ages) * 8))

names(opt.both.inctest) <- c("Incidence", "Age", "Group")
names(opt.both.prevtest) <- c("Prevalence_uncond", "Prevalence_cond", "Age", "State", "Source")






# matrices for incidence and prevalence optimizations
# opt.prev.matrices <- make_trans_matrix(opt.prev$solution, r45.params)
# opt.inc.matrices <- make_trans_matrix(opt.inc$solution, r45.params)
# 
# 
# # show how optimization with prevalence does well with prev, bad with incidence
# opt.prev.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014, 1.65, 1, 
#                                         k0 = opt.prev.matrices[[1]], k1 = opt.prev.matrices[[2]]) * 100
# opt.prev.prevalence <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1, 
#                                                    k0 = opt.prev.matrices[[1]], opt.prev.matrices[[2]])
# 
# # show how optimization with incidence does well with incidence, bad with prev
# opt.inc.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014, 1.65, 1, 
#                                        k0 = opt.inc.matrices[[1]], k1 = opt.inc.matrices[[2]]) * 100
# opt.inc.prevalence <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
#                                                   k0 = opt.inc.matrices[[1]], k1 = opt.inc.matrices[[2]])
# 
# inc.imperfect <- cbind.data.frame(c(opt.prev.incidence, opt.inc.incidence, empirical.incidence),
#                               rep(inc.ages, 3),
#                               rep(c("Prevalence-Optimized", "Incidence-Optimized", "Jack"), each = length(inc.ages)))
# 
# prev.imperfect <- cbind.data.frame(c(as.vector(opt.prev.prevalence), as.vector(opt.inc.prevalence), as.vector(fem_prev_u)),
#                                    rep(ages, 3 * 8),
#                                    rep(c("Prevalence-Optimized", "Incidence-Optimized", "Jack"), each = 8 * length(ages)),
#                                    rep(rep(c("State 1", "State 2", "State 3", "State 4",
#                                              "State 6", "State 7", "State 8", "State 9"), each = length(ages)), 3))
# 
# names(inc.imperfect) <- c("Incidence", "Age", "Group")
# names(prev.imperfect) <- c("Prevalence", "Age", "Group", "State")
# 
# ggplot(inc.imperfect, mapping = aes(x = Age, y = Incidence, color = Group)) +
#   geom_line(size = .75) +
#   scale_y_continuous(trans = "log2", breaks = c(.05, .1, .25, .5, 1, 2, 5, 10)) + 
#   scale_color_manual(values = c("black", "green", "red")) + 
#   labs(title = "Empirical vs Multistate Incidence of AD Dementia",
#        y = "Incidence (% per year)")
# ggsave("imperfect.incidence.pdf")
# ggplot(prev.imperfect, mapping = aes(x = Age, y = Prevalence, color = Group)) +
#   geom_line(size = .75) + 
#   facet_wrap( ~ State) + 
#   labs(title = "Jack vs Multistate unconditional prevalence estimates",
#        y = "Unconditional Prevalence (%)")
# ggsave("imperfect.prevalence.pdf")
# 
# # saving parameter estimates and data for prev/inc comparison
# save(opt.both.prevtest, file = "dual.opt.prevalence.R")
# save(opt.both.inctest, file = "dual.opt.incidence.R")
# 
# opt.both.params <- as.data.frame(opt.both$solution)
# save(opt.both.params, file = "dual.opt.parameters.R")
# save(inc.imperfect, file = "inc.imperfect.R")
# save(prev.imperfect, file = "prev.imperfect.R")

# it hadn't converged before with 15k iterations, so let's run again til convergence
# convergence is in terms of relative tolerance of x

opt.both.rerun <- nloptr(x0 = opt.both$solution,
                         eval_f = eval_f_both,
                         lb = lb, ub = ub,
                         eval_g_ineq = eval_g_ineq_both,
                         opts = list("algorithm"="NLOPT_LN_COBYLA",
                                     "xtol_rel"=1e-3,
                                     "maxeval"=10000),
                         r45 = r45.params,
                         prevs = fem_prev_u,
                         incidence = empirical.incidence)


opt.both.params_03.19 <- as.data.frame(opt.both.rerun$solution)
save(opt.both.params_03.19, file = "dual.opt.parameters_03.19.2021.R")

opt.both.matrices <- make_trans_matrix(opt.both.rerun$solution, r45.params)

# incidence, unconditional prevalence, and conditional prevalence estimates using optimized parameters
opt.both.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014,
                                        k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]]) * 100
opt.both.prevalence.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
                                                     k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])
opt.both.prevalence.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1,
                                              k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])

# datasets for comparing with empirical incidence, jack prevalence
opt.both.inctest <- cbind.data.frame(c(opt.both.incidence, empirical.incidence),
                                     rep(inc.ages, 2),
                                     c(rep("Dual-optimized", length(inc.ages)),
                                       rep("Empirical", length(inc.ages))))
opt.both.prevtest <- cbind.data.frame(c(as.vector(opt.both.prevalence.u), as.vector(fem_prev_u)),
                                      c(as.vector(opt.both.prevalence.c), as.vector(fem_prev)),
                                      rep(ages, 16),
                                      rep(rep(c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                "T+", "T+ N+", "N+", "A+ N+"), each = length(ages)), 2),
                                      rep(c("Dual-optimized", "Jack"), each = length(ages) * 8))

names(opt.both.inctest) <- c("Incidence", "Age", "Group")
names(opt.both.prevtest) <- c("Prevalence_uncond", "Prevalence_cond", "Age", "State", "Source")

opt.both.prevtest$State <- factor(opt.both.prevtest$State, levels = c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                                      "T+", "T+ N+", "N+", "A+ N+"))

save(opt.both.inctest, file = "dual.opt.incidence_03.19.2021.R")
save(opt.both.prevtest, file = "dual.opt.prevalence_03.19.2021.R")


base.loss.both <- eval_f_both_logs(opt.both.rerun$solution, r45.params, fem_prev_u, empirical.incidence)
base.logs <- eval_f_both_logs(init.simple, r45.params, fem_prev_u, empirical.incidence)
base.nologs <- eval_f_both_nologs(init.simple, r45.params, fem_prev_u, empirical.incidence)
opt.both.logs <- nloptr(x0 = init.simple, 
                        eval_f = eval_f_both_logs, 
                        lb = lb, ub = ub, 
                        eval_g_ineq = eval_g_ineq_both,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1e-3,
                                    "maxeval"=30000),
                        r45 = r45.params,
                        prevs = fem_prev_u,
                        incidence = empirical.incidence)




########### Optimization giving equal weight to incidence and prevalence losses

ew.base.1 <- eval_f_both_logs(init.simple, r45.params, fem_prev_u, empirical.incidence)
ew.base.2 <- eval_f_both_logs(init.simple2, r45.params, fem_prev_u, empirical.incidence)
opt.equal.weight <- nloptr(x0 = init.simple, 
                        eval_f = eval_f_both_logs, 
                        lb = lb, ub = ub, 
                        eval_g_ineq = eval_g_ineq_both,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1e-3,
                                    "maxeval"=25000),
                        r45 = r45.params,
                        prevs = fem_prev_u,
                        incidence = empirical.incidence)

init.simple2 <- rep(c(-9, 0.05), 12)
opt.ew.2 <- nloptr(x0 = init.simple2 , 
                   eval_f = eval_f_both_logs, 
                   lb = lb, ub = ub, 
                   eval_g_ineq = eval_g_ineq_both,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=25000),
                   r45 = r45.params,
                   prevs = fem_prev_u,
                   incidence = empirical.incidence)



opt.ew1.matrices <- make_trans_matrix(opt.equal.weight$solution, r45.params)
opt.ew2.matrices <- make_trans_matrix(opt.ew.2$solution, r45.params)

ew1.inc <- incidence.f.multi(inc.ages, 2014, intalpha, 2014,
                             k0 = opt.ew1.matrices[[1]], k1 = opt.ew1.matrices[[2]]) * 100
ew2.inc <- incidence.f.multi(inc.ages, 2014, intalpha, 2014,
                             k0 = opt.ew2.matrices[[1]], k1 = opt.ew2.matrices[[2]]) * 100
ew1.prev.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
                                                     k0 = opt.ew1.matrices[[1]], k1 = opt.ew1.matrices[[2]])
ew2.prev.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
                                          k0 = opt.ew2.matrices[[1]], k1 = opt.ew2.matrices[[2]])
ew1.prev.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1,
                                          k0 = opt.ew1.matrices[[1]], k1 = opt.ew1.matrices[[2]])
ew2.prev.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1,
                                          k0 = opt.ew2.matrices[[1]], k1 = opt.ew2.matrices[[2]])

ew.inctest <- cbind.data.frame(c(ew1.inc, empirical.incidence),
                                     rep(inc.ages, 2),
                                     c(rep("Inits 1", length(inc.ages)),
                                       rep("Empirical", length(inc.ages))))
ew.prevtest <- cbind.data.frame(c(as.vector(ew1.prev.u), as.vector(fem_prev_u)),
                                      c(as.vector(ew1.prev.c), as.vector(fem_prev)),
                                      rep(ages, 8 * 2),
                                      rep(rep(c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                "T+", "T+ N+", "N+", "A+ N+"), each = length(ages)), 2),
                                      rep(c("Inits 1", "Jack"), each = length(ages) * 8))

names(ew.inctest) <- c("Incidence", "Age", "Group")
names(ew.prevtest) <- c("Prevalence_uncond", "Prevalence_cond", "Age", "State", "Source")

ew.prevtest$State <- factor(ew.prevtest$State, levels = c("Normal", "A+", "A+ T+", "A+ T+ N+",
                                                                      "T+", "T+ N+", "N+", "A+ N+"))


plot.uncond <- ggplot(data = ew.prevtest, mapping = aes(x = Age, y = Prevalence_uncond)) + 
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Unconditional Prevalence",
       y = "Unconditional Prevalence (%)")

plot.cond <- ggplot(data = ew.prevtest, mapping = aes(x = Age, y = Prevalence_cond)) +
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Conditional Prevalence",
       y = "Conditional Prevalence (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25))

plot.inc <- ggplot(ew.inctest, mapping = aes(x = Age, y = Incidence, color = Group)) +
  geom_line(size = .75) +
  scale_y_continuous(trans = "log2", breaks = c(.05, .1, .25, .5, 1, 2, 5, 10)) + 
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Empirical and Fitted Incidence of AD Dementia",
       y = "Incidence (% per year)")






# incidence, unconditional prevalence, and conditional prevalence estimates using optimized parameters
opt.both.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014,
                                        k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]]) * 100
opt.both.prevalence.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1,
                                                     k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])
opt.both.prevalence.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1,
                                              k0 = opt.both.matrices[[1]], k1 = opt.both.matrices[[2]])

# datasets for comparing with empirical incidence, jack prevalence
opt.both.inctest <- cbind.data.frame(c(opt.both.incidence, empirical.incidence),
                                     rep(inc.ages, 2),
                                     c(rep("Dual-optimized", length(inc.ages)),
                                       rep("Empirical", length(inc.ages))))
opt.both.prevtest <- cbind.data.frame(c(as.vector(opt.both.prevalence.u), as.vector(fem_prev_u)),
                                      c(as.vector(opt.both.prevalence.c), as.vector(fem_prev)),
                                      rep(ages, 16),
                                      rep(rep(c("State 1", "State 2", "State 3", "State 4",
                                                "State 6", "State 7", "State 8", "State 9"), each = length(ages)), 2),
                                      rep(c("Dual-optimized", "Jack"), each = length(ages) * 8))



opt.logs.matrices <- make_trans_matrix(opt.both.logs$solution, r45.params)



#### Results from dual optimization with both log(prev) and log(inc)

# opt.logs.incidence <- incidence.f.multi(inc.ages, 2014, intalpha, 2014, 
#                                         k0 = opt.logs.matrices[[1]], k1 = opt.logs.matrices[[2]]) * 100
# opt.logs.prevalence.u <- Prevrate.f.multi.ATN.uncond(ages, 2014, intalpha, 2014, 1.65, 1, 
#                                                      k0 = opt.logs.matrices[[1]], k1 = opt.logs.matrices[[2]])
# opt.logs.prevalence.c <- Prevrate.f.multi.ATN(ages, 2014, intalpha, 2014, 1.65, 1, 
#                                               k0 = opt.logs.matrices[[1]], k1 = opt.logs.matrices[[2]])
# 
# 
# opt.logs.inctest <- cbind.data.frame(c(opt.logs.incidence, opt.both.incidence, empirical.incidence),
#                                        rep(inc.ages, 3),
#                                        c(rep("Logged-both", length(inc.ages)),
#                                          rep("Logged-incidence", length(inc.ages)),
#                                          rep("Empirical", length(inc.ages))))
# opt.logs.prevtest <- cbind.data.frame(c(as.vector(opt.logs.prevalence.u), as.vector(opt.both.prevalence.u), as.vector(fem_prev_u)),
#                                         c(as.vector(opt.logs.prevalence.c), as.vector(opt.both.prevalence.c), as.vector(fem_prev)),
#                                         rep(ages, 3 * 8),
#                                         rep(rep(c("State 1", "State 2", "State 3", "State 4",
#                                                   "State 6", "State 7", "State 8", "State 9"), each = length(ages)), 3),
#                                         rep(c("Logged-both", "Logged-incidence", "Jack"), each = length(ages) * 8))
# names(opt.logs.inctest) <- c("Incidence", "Age", "Group")
# names(opt.logs.prevtest) <- c("Prevalence_uncond", "Prevalence_cond", "Age", "State", "Source")
# 
# save(opt.logs.inctest, file = "opt.logs.inc_04.08.2021.R")
# save(opt.logs.prevtest, file = "opt.logs.prev_04.08.2021.R")
# save(opt.logs.params, file = "opt.logs.params_04.08.2021.R")

base.weighted.2 <- eval_f_logs_weighted(init.simple, r45 = r45.params, prevs = fem_prev_u, incidence = empirical.incidence, w = 1)

inits.simple <- matrix(c(init.simple, rep(c(-6, 0.04), 12)), byrow = FALSE, ncol = 2)
weighted.opts <- list()

for(i in 1:2){
  weighted.opts[[i]] <- nloptr(x0 = inits.simple[,i], 
                          eval_f = eval_f_logs_weighted, 
                          lb = lb, ub = ub, 
                          eval_g_ineq = eval_g_ineq_weighted,
                          opts = list("algorithm"="NLOPT_LN_COBYLA",
                                      "xtol_rel"=1e-3,
                                      "maxeval"=10000),
                          r45 = r45.params,
                          prevs = fem_prev_u,
                          incidence = empirical.incidence,
                          w = 2)
}
inc.weighted.21 <- nloptr(x0 = inits.simple[,2], 
                         eval_f = eval_f_logs_weighted, 
                         lb = lb, ub = ub, 
                         eval_g_ineq = eval_g_ineq_weighted,
                         opts = list("algorithm"="NLOPT_LN_COBYLA",
                                     "xtol_rel"=1e-3,
                                     "maxeval"=10),
                         r45 = r45.params,
                         prevs = fem_prev_u,
                         incidence = empirical.incidence,
                         w = 2)




weight.2.prev <- make_prevplot_data(inc.weighted.2$solution)
weight.2.inc <- make_incplot_data(inc.weighted.2$solution)
weight.2.params <- matrix(inc.weighted.2$solution, ncol = 2, byrow = TRUE)

save(weight.2.prev, file = "weight.2.prevs_04.21.2021.R")
save(weight.2.inc, file = "weight.2.inc_04.21.2021.R")
saveRDS(inc.weighted.2, file = "weight.2.opt_04.21.2021.rds")
saveRDS(opt.equal.weight, file = "weight.1.opt_04.21.2021.rds")

weight.2.uncond.plot <- ggplot(data = weight.2.prev, mapping = aes(x = Age, y = Prevalence_uncond)) + 
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Unconditional Prevalence",
       y = "Unconditional Prevalence (%)")

weight.2.cond.plot <- ggplot(data = weight.2.prev, mapping = aes(x = Age, y = Prevalence_cond)) +
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Conditional Prevalence",
       y = "Conditional Prevalence (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25))

weight.2.inc.plot <- ggplot(weight.2.inc, mapping = aes(x = Age, y = Incidence, color = Source)) +
  geom_line(size = .75) +
  scale_y_continuous(trans = "log2", breaks = c(.05, .1, .25, .5, 1, 2, 5, 10)) + 
  scale_color_manual(values = c("black", "deepskyblue1")) +
  labs(title = "Empirical and Fitted Incidence of AD Dementia",
       y = "Incidence (% per year)")


##### optimize using averaged male/female prevalence

avg.prev.opt <- nloptr(x0 = init.simple, 
                          eval_f = eval_f_logs_weighted, 
                          lb = lb, ub = ub, 
                          eval_g_ineq = eval_g_ineq_weighted,
                          opts = list("algorithm"="NLOPT_LN_COBYLA",
                                      "xtol_rel"=1e-3,
                                      "maxeval"=12000),
                          r45 = r45.params,
                          prevs = avg_prev_u,
                          incidence = empirical.incidence,
                          w = 1)

saveRDS(avg.prev.opt, file = "avg.prev.opt_05.03.2021.rds")
avg.prev <- make_prevplot_data(avg.prev.opt$solution)
avg.inc <- make_incplot_data(avg.prev.opt$solution)

avg.uncond.plot <- ggplot(data = avg.prev, mapping = aes(x = Age, y = Prevalence_uncond)) + 
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Unconditional Prevalence",
       y = "Unconditional Prevalence (%)")

avg.cond.plot <- ggplot(data = avg.prev, mapping = aes(x = Age, y = Prevalence_cond)) +
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Conditional Prevalence",
       y = "Conditional Prevalence (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25))

avg.inc.plot <- ggplot(avg.inc, mapping = aes(x = Age, y = Incidence, color = Source)) +
  geom_line(size = .75) +
  scale_color_manual(values = c("black", "deepskyblue1")) +
  labs(title = "Empirical and Fitted Incidence of AD Dementia",
       y = "Incidence (% per year)")


#### tweaking the optimal solution a little so that when it searches
#### it doesn't go outside the range of constraints

ap.opt.rerun <- nloptr(x0 = init.simple, 
                       eval_f = eval_f_logs_weighted, 
                       lb = lb, ub = ub, 
                       eval_g_ineq = eval_g_ineq_weighted,
                       opts = list("algorithm"="NLOPT_LN_COBYLA",
                                   "xtol_rel"=1e-3,
                                   "ftol_rel"=1e-5,
                                   "maxeval"=20000),
                       r45 = r45.params,
                       prevs = avg_prev_u,
                       incidence = empirical.incidence,
                       w = 1)

saveRDS(ap.opt.rerun, file = "avg.prev.opt_05.05.2021.rds")
avg.prev <- make_prevplot_data(avg.prev.opt$solution)
avg.inc <- make_incplot_data(avg.prev.opt$solution)

avg.uncond.plot <- ggplot(data = avg.prev, mapping = aes(x = Age, y = Prevalence_uncond)) + 
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Unconditional Prevalence",
       y = "Unconditional Prevalence (%)")

avg.cond.plot <- ggplot(data = avg.prev, mapping = aes(x = Age, y = Prevalence_cond)) +
  geom_line(aes(color = Source)) + 
  facet_wrap( ~ State, nrow = 4) +
  scale_color_manual(values = c("deepskyblue1", "black")) +
  labs(title = "Cross-sectional and Fitted Conditional Prevalence",
       y = "Conditional Prevalence (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25))

avg.inc.plot <- ggplot(avg.inc, mapping = aes(x = Age, y = Incidence, color = Source)) +
  geom_line(size = .75) +
  scale_color_manual(values = c("black", "deepskyblue1")) +
  labs(title = "Empirical and Fitted Incidence of AD Dementia",
       y = "Incidence (% per year)")


####### Seeing how much estimates change with more iterations
opt.500 <- nloptr(x0 = init.simple, 
                   eval_f = eval_f_logs_weighted, 
                   lb = lb, ub = ub, 
                   eval_g_ineq = eval_g_ineq_weighted,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=500),
                   r45 = r45.params,
                   prevs = avg_prev_u,
                   incidence = empirical.incidence,
                   w = 1)

opt.1000 <- nloptr(x0 = init.simple, 
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
opt.3000 <- nloptr(x0 = init.simple, 
                  eval_f = eval_f_logs_weighted, 
                  lb = lb, ub = ub, 
                  eval_g_ineq = eval_g_ineq_weighted,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1e-3,
                              "maxeval"=3000),
                  r45 = r45.params,
                  prevs = avg_prev_u,
                  incidence = empirical.incidence,
                  w = 1)
opt.5000 <- nloptr(x0 = init.simple, 
                   eval_f = eval_f_logs_weighted, 
                   lb = lb, ub = ub, 
                   eval_g_ineq = eval_g_ineq_weighted,
                   opts = list("algorithm"="NLOPT_LN_COBYLA",
                               "xtol_rel"=1e-3,
                               "maxeval"=5000),
                   r45 = r45.params,
                   prevs = avg_prev_u,
                   incidence = empirical.incidence,
                   w = 1)

# write down these results

saveRDS(opt.500, file = "OptResults/simple.500.rds")
saveRDS(opt.1000, file = "OptResults/simple.1000.rds")
saveRDS(opt.3000, file = "OptResults/simple.3000.rds")
saveRDS(opt.5000, file = "OptResults/simple.5000.rds")

