# Assessing the potential effects of drug Aducanumab on
# lifetime/10-year risks and prevalence of AD dementia

library(tidyverse)
library(here)

source(here("BrookFuncs.R"))

drug.effects <- c(1, 0.95, 0.9, 0.75)

alpha.drug.effects <- list()
for(i in 1:length(drug.effects)){
  alpha.drug.effects[[i]] <- intalpha.AN
  alpha.drug.effects[[i]][5, 7] <- drug.effects[i]
}

# list of lifetime tables for Aducanumab relative risks

Aducan.lifetime.f <- list()
Aducan.lifetime.m <- list()
Aducan.tenyear.f <- list()
Aducan.tenyear.m <- list()

lifetime.ages <- seq(60, 95, 5)

delay <- 2

for(i in 1:length(alpha.drug.effects)){
  Aducan.lifetime.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.lifetime.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.tenyear.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.tenyear.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  
  for(j in 1:length(lifetime.ages)){
    for(k in 1:6){
      lifetime.f.curr <- lifetime.AN(lifetime.ages[j], "Female", k, 2017 + delay, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      lifetime.m.curr <- lifetime.AN(lifetime.ages[j], "Male", k, 2017 + delay, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      Aducan.lifetime.f[[i]][j, k] <- lifetime.f.curr[1]
      Aducan.lifetime.m[[i]][j, k] <- lifetime.m.curr[1]
      Aducan.tenyear.f[[i]][j, k] <-  lifetime.f.curr[2]
      Aducan.tenyear.m[[i]][j, k] <-  lifetime.m.curr[2]
      
      Aducan.lifetime.f[[i]] <-  Aducan.lifetime.f[[i]]
      Aducan.lifetime.m[[i]] <-  Aducan.lifetime.m[[i]]
      Aducan.tenyear.f[[i]] <- Aducan.tenyear.f[[i]]
      Aducan.tenyear.m[[i]] <- Aducan.tenyear.m[[i]]
      
    }
  }
}

for(i in 1:length(Aducan.lifetime.f)){
  Aducan.lifetime.f[[i]] <- cbind(seq(60, 95, 5), Aducan.lifetime.f[[i]])
  Aducan.lifetime.m[[i]] <- cbind(seq(60, 95, 5), Aducan.lifetime.m[[i]])
  Aducan.tenyear.f[[i]] <- cbind(seq(60, 95, 5), Aducan.tenyear.f[[i]])
  Aducan.tenyear.m[[i]] <- cbind(seq(60, 95, 5), Aducan.tenyear.m[[i]])
  
  colnames(Aducan.lifetime.f[[i]]) <- colnames(Aducan.lifetime.m[[i]]) <- 
    colnames(Aducan.tenyear.f[[i]]) <- colnames(Aducan.tenyear.m[[i]]) <- 
    c("Age", "Normal[1]", "A[2]", "N[3]", "A+N[4]", "A+N+MCI[5]", "N+MCI[6]")
  rownames(Aducan.lifetime.f[[i]]) <- rownames(Aducan.lifetime.m[[i]]) <- 
    rownames(Aducan.tenyear.f[[i]]) <- rownames(Aducan.tenyear.m[[i]]) <- 
    NULL
}

names(Aducan.lifetime.f) <- names(Aducan.lifetime.m) <- names(Aducan.tenyear.f) <- names(Aducan.tenyear.m) <-
  c("RR=1", "RR=0.95", "RR=0.9", "RR=0.75")

Aducan.all <- list(Aducan.lifetime.f, Aducan.lifetime.m, Aducan.tenyear.f, Aducan.tenyear.m)

names(Aducan.all) <- c("lifetime.f", "lifetime.m", "tenyear.f", "tenyear.m")

# saveRDS(Aducan.all, "Aducanumab/Aducan.lifetime.rds")

Aducan.nodelay.lifetime.f <- list()
Aducan.nodelay.lifetime.m <- list()
Aducan.nodelay.tenyear.f <- list()
Aducan.nodelay.tenyear.m <- list()

lifetime.ages <- seq(60, 95, 5)

for(i in 1:length(alpha.drug.effects)){
  Aducan.nodelay.lifetime.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.nodelay.lifetime.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.nodelay.tenyear.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.nodelay.tenyear.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  
  for(j in 1:length(lifetime.ages)){
    for(k in 1:6){
      lifetime.f.curr <- lifetime.AN(lifetime.ages[j], "Female", k, 2017, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      lifetime.m.curr <- lifetime.AN(lifetime.ages[j], "Male", k, 2017, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      Aducan.nodelay.lifetime.f[[i]][j, k] <- cbind(seq(60, 95, 5), lifetime.f.curr[1])
      Aducan.nodelay.lifetime.m[[i]][j, k] <- cbind(seq(60, 95, 5), lifetime.m.curr[1])
      Aducan.nodelay.tenyear.f[[i]][j, k] <- cbind(seq(60, 95, 5), lifetime.f.curr[2])
      Aducan.nodelay.tenyear.m[[i]][j, k] <- cbind(seq(60, 95, 5), lifetime.m.curr[2])
      
      colnames(Aducan.nodelay.lifetime.f[[i]]) <- colnames(Aducan.nodelay.lifetime.m[[i]]) <- 
        colnames(Aducan.nodelay.tenyear.f[[i]]) <- colnames(Aducan.nodelay.tenyear.m[[i]]) <- 
        c("Age", "Normal[1]", "A[2]", "N[3]", "A+N[4]", "A+N+MCI[5]", "N+MCI[6]")
      
      rownames(Aducan.nodelay.lifetime.f[[i]]) <- rownames(Aducan.nodelay.lifetime.m[[i]]) <- 
        rownames(Aducan.nodelay.tenyear.f[[i]]) <- rownames(Aducan.nodelay.tenyear.m[[i]]) <- 
        NULL
      
    }
  }
}


names(Aducan.nodelay.lifetime.f) <- names(Aducan.nodelay.lifetime.m) <- names(Aducan.nodelay.tenyear.f) <- names(Aducan.nodelay.tenyear.m) <-
  c("RR=1", "RR=0.95", "RR=0.9", "RR=0.75")

Aducan.nodelay.all <- list(Aducan.nodelay.lifetime.f, Aducan.nodelay.lifetime.m, Aducan.nodelay.tenyear.f, Aducan.nodelay.tenyear.m)

names(Aducan.nodelay.all) <- c("lifetime.f", "lifetime.m", "tenyear.f", "tenyear.m")

# saveRDS(Aducan.nodelay.all, "Aducanumab/Aducan.nodelay.lifetime.rds")

Aducan.lifetime.f <- list()
Aducan.lifetime.m <- list()
Aducan.tenyear.f <- list()
Aducan.tenyear.m <- list()

lifetime.ages <- seq(60, 95, 5)

delay <- 2

for(i in 1:length(alpha.drug.effects)){
  Aducan.lifetime.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.lifetime.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.tenyear.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  Aducan.tenyear.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 6)
  
  for(j in 1:length(lifetime.ages)){
    for(k in 1:6){
      lifetime.f.curr <- lifetime.AN(lifetime.ages[j], "Female", k, 2017 + delay, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      lifetime.m.curr <- lifetime.AN(lifetime.ages[j], "Male", k, 2017 + delay, alpha = alpha.drug.effects[[i]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      Aducan.lifetime.f[[i]][j, k] <- lifetime.f.curr[1]
      Aducan.lifetime.m[[i]][j, k] <- lifetime.m.curr[1]
      Aducan.tenyear.f[[i]][j, k] <-  lifetime.f.curr[2]
      Aducan.tenyear.m[[i]][j, k] <-  lifetime.m.curr[2]
      
      Aducan.lifetime.f[[i]] <-  Aducan.lifetime.f[[i]]
      Aducan.lifetime.m[[i]] <-  Aducan.lifetime.m[[i]]
      Aducan.tenyear.f[[i]] <- Aducan.tenyear.f[[i]]
      Aducan.tenyear.m[[i]] <- Aducan.tenyear.m[[i]]
      
    }
  }
}

for(i in 1:length(Aducan.lifetime.f)){
  Aducan.lifetime.f[[i]] <- cbind(seq(60, 95, 5), Aducan.lifetime.f[[i]])
  Aducan.lifetime.m[[i]] <- cbind(seq(60, 95, 5), Aducan.lifetime.m[[i]])
  Aducan.tenyear.f[[i]] <- cbind(seq(60, 95, 5), Aducan.tenyear.f[[i]])
  Aducan.tenyear.m[[i]] <- cbind(seq(60, 95, 5), Aducan.tenyear.m[[i]])
  
  colnames(Aducan.lifetime.f[[i]]) <- colnames(Aducan.lifetime.m[[i]]) <- 
    colnames(Aducan.tenyear.f[[i]]) <- colnames(Aducan.tenyear.m[[i]]) <- 
    c("Age", "Normal[1]", "A[2]", "N[3]", "A+N[4]", "A+N+MCI[5]", "N+MCI[6]")
  rownames(Aducan.lifetime.f[[i]]) <- rownames(Aducan.lifetime.m[[i]]) <- 
    rownames(Aducan.tenyear.f[[i]]) <- rownames(Aducan.tenyear.m[[i]]) <- 
    NULL
}

names(Aducan.lifetime.f) <- names(Aducan.lifetime.m) <- names(Aducan.tenyear.f) <- names(Aducan.tenyear.m) <-
  c("RR=1", "RR=0.95", "RR=0.9", "RR=0.75")

Aducan.all <- list(Aducan.lifetime.f, Aducan.lifetime.m, Aducan.tenyear.f, Aducan.tenyear.m)

names(Aducan.all) <- c("lifetime.f", "lifetime.m", "tenyear.f", "tenyear.m")

# saveRDS(Aducan.all, "Aducanumab/Aducan.lifetime.rds")

########## One more analysis but with big effects

# matrices will have a different format
# one row per age, but matrix will be for a state instead of a RR

Aducan.big.lifetime.f <- list()
Aducan.big.lifetime.m <- list()
Aducan.big.tenyear.f <- list()
Aducan.big.tenyear.m <- list()

big.effects <- c(1, 0.75, 0.5, 0.25)

alpha.big.effects <- list()
for(i in 1:length(drug.effects)){
  alpha.big.effects[[i]] <- intalpha.AN
  alpha.big.effects[[i]][5, 7] <- big.effects[i]
}

lifetime.ages <- seq(60, 95, 5)

for(i in 1:5){
  Aducan.big.lifetime.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 5)
  Aducan.big.lifetime.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 5)
  Aducan.big.tenyear.f[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 5)
  Aducan.big.tenyear.m[[i]] <- matrix(nrow = length(lifetime.ages), ncol = 5)
  
  for(j in 1:length(lifetime.ages)){
    for(k in 1:length(big.effects)){
      lifetime.f.curr <- lifetime.AN(lifetime.ages[j], "Female", i, 2017, alpha = alpha.big.effects[[k]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      lifetime.m.curr <- lifetime.AN(lifetime.ages[j], "Male", i, 2017, alpha = alpha.big.effects[[k]], 
                                     k0 = k0.AN, k1 = k1.AN, k012 = k012.AN, k112 = k112.AN, k12 = k12.AN)
      Aducan.big.lifetime.f[[i]][j, k + 1] <- lifetime.f.curr[1]
      Aducan.big.lifetime.m[[i]][j, k + 1] <- lifetime.m.curr[1]
      Aducan.big.tenyear.f[[i]][j, k + 1] <- lifetime.f.curr[2]
      Aducan.big.tenyear.m[[i]][j, k + 1] <- lifetime.m.curr[2]
      
    }
  }
  Aducan.big.lifetime.f[[i]][,1] <- Aducan.big.lifetime.m[[i]][,1] <- 
    Aducan.big.tenyear.f[[i]][,1] <- Aducan.big.tenyear.m[[i]][,1] <- seq(60, 95, 5)
  
  colnames(Aducan.big.lifetime.f[[i]]) <- colnames(Aducan.big.lifetime.m[[i]]) <- 
    colnames(Aducan.big.tenyear.f[[i]]) <- colnames(Aducan.big.tenyear.m[[i]]) <- 
    c("Age", "RR=1", "RR=0.75", "RR=0.5", "RR=0.25")
  rownames(Aducan.big.lifetime.f[[i]]) <- rownames(Aducan.big.lifetime.m[[i]]) <- 
    rownames(Aducan.big.tenyear.f[[i]]) <- rownames(Aducan.big.tenyear.m[[i]]) <- 
    NULL
}

names(Aducan.big.lifetime.f) <- names(Aducan.big.lifetime.m) <- 
  names(Aducan.big.tenyear.f) <- names(Aducan.big.tenyear.m) <- c("state1", "state2", "state3", "state4", "state5")

big.all <- list(Aducan.big.lifetime.f, Aducan.big.lifetime.m,
                Aducan.big.tenyear.f, Aducan.big.tenyear.m)
names(big.all) <- c("lifetime.f", "lifetime.m", "tenyear.f", "tenyear.m")

saveRDS(big.all, file = "Aducanumab/big.effects.rds")
