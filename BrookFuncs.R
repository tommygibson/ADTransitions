####### This will hold all the functions being used in the alzheimer's research

### First section is ATN model with no piecewise component
### Second component starts ~ line 260 and is AN model (Nada's old code)

library(readxl)
library(here)


##### data setup, known stuff
##### includes death rates, prevalence rates from jack (2017), incidence rates from meta-analysis
# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Brookmeyer/Code")
rates <- read.csv(here("U.S.mortality.rates.csv"))
jack_prev_men <- as.matrix(read_excel(here("ATN_estimates_rescaled.xlsx"))[,c(1, 2, 6, 8, 9, 4, 5, 3, 7)])
men_prev <- jack_prev_men[,-1] / 100
jack_prev_fem <- as.matrix(read_excel(here("ATN_estimates_rescaled.xlsx"))[,c(1, 10, 14, 16, 17, 12, 13, 11, 15)])
fem_prev <- jack_prev_fem[,-1] / 100
avg_prev <- (fem_prev + men_prev) / 2

prev.ages <- 50:95
inc.ages <- 65:90

empirical.incidence <- 0.00117 * exp(0.126 * (inc.ages - 60)) * 100


# trans.params <- read.csv("params.opt_01.26.2020.csv")[,-1]
r45.params <- unlist(read.csv(here("params.r45.csv"))[-1])
r45.params[1] <- exp(r45.params[1])

# p.preclinical <- as.vector(read.csv("prev.preclinical_03.01.2021.csv")[,2])
# p.preclinical <- as.vector(read.csv(here("prev.preclinical_03.15.2021.csv"))[,2])
p.preclinical <- as.vector(read.csv(here("prev.preclinical_06.10.2021.csv"))[,2])

fem_prev_u <- p.preclinical * fem_prev
men_prev_u <- p.preclinical * men_prev

avg_prev_u <- (fem_prev_u + men_prev_u) / 2

death <- 12
ad.state <- 10

#females
dr.f <- function(a, y, f){
  rates.sub <- subset(rates, rates[, 1] == y)
  d.f <- f * rates.sub[a + 1, 3]
  return(d.f)
}
#males
dr.m <- function(a, y, f){
  rates.sub <- subset(rates, rates[, 1] == y)
  d.f <- f * rates.sub[a + 1, 4]
  return(d.f)
}

alpha.int<-function(one.two, one.six, one.eight, two.three, two.nine, three.four, four.five, five.ten, six.three,
                    six.seven, seven.four, eight.seven, eight.nine, nine.four, ten.eleven){
  
  # 12 states including advanced alz (11) and death (12)
  alpha<-matrix(1, nrow = 11, ncol = 12)
  
  alpha[1, 2] <- one.two
  alpha[1, 6] <- one.six
  alpha[1, 8] <- one.eight
  alpha[2, 3] <- two.three
  alpha[2, 9] <- two.nine
  alpha[3, 4] <- three.four
  alpha[4, 5] <- four.five
  alpha[5, 10] <- five.ten
  alpha[6, 3] <- six.three
  alpha[6, 7] <- six.seven
  alpha[7, 4] <- seven.four
  alpha[8, 7] <- eight.seven
  alpha[8, 9] <- eight.nine
  alpha[9, 4] <- nine.four
  alpha[10, 11] <- ten.eleven  ## stage 11 is advanced alzheimer's
  return(alpha)
}

intalpha = do.call(alpha.int, as.list(rep(1, 15)))

alpha.int.reduced <-function(one.two, one.six, two.three, two.seven, three.four, four.five, five.eight, 
                             six.three, six.four, six.seven, seven.four, eight.nine){
  
  # 12 states including advanced alz (9) and death (10)
  alpha <- matrix(1, nrow = 9, ncol = 10)
  
  alpha[1, 2] <- one.two
  alpha[1, 6] <- one.six
  alpha[2, 3] <- two.three
  alpha[2, 7] <- two.seven
  alpha[3, 4] <- three.four
  alpha[4, 5] <- four.five
  alpha[5, 8] <- five.eight
  alpha[6, 3] <- six.three
  alpha[6, 4] <- six.four
  alpha[6, 7] <- six.seven
  alpha[7, 4] <- seven.four
  alpha[8, 9] <- eight.nine  ## stage 9 is advanced alzheimer's
  return(alpha)
}

intalpha.reduced <- do.call(alpha.int.reduced, as.list(rep(1, 12)))


alpha.nine<-function(one.two,  two.four, four.five, one.three, three.four,three.six, 
                     five.seven, six.seven, seven.eight){
  alpha<-matrix(1,nrow=8,ncol=9)
  alpha[1,2]<-one.two
  alpha[2,4]<-two.four
  alpha[4,5]<-four.five
  alpha[1,3]<-one.three
  alpha[3,4]<-three.four
  alpha[3,6]<-three.six
  alpha[5,7]<-five.seven
  alpha[6,7]<-six.seven
  alpha[7,8]<-seven.eight
  return(alpha)
}
#example: no intervention 
intalpha.AN = alpha.nine(1,1,1,1,1,1,1,1,1)

########### FUNCTION MAKES MATRICES k0 and k1 FROM A VECTOR OF (log(k0ij), k1ij) pairs
# vector length is 24
make_trans_matrix <- function(x, r45){
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
  
  return(list(k0, k1))
}

make_trans_matrix_high <- function(x, r45){
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
  
  return(list(k0, k1))
}

make_trans_matrix_low <- function(x, r45){
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
  
  return(list(k0, k1))
}

make_trans_matrix_reduced <- function(x, r45){
  k0 <- matrix(0, nrow = 7, ncol = 8)
  k1 <- matrix(0, nrow = 7, ncol = 8)
  
  x[seq(1, 17, 2)] <- exp(x[seq(1, 17, 2)])
  
  # from state 1
  k0[1, 2] <- x[1]
  k1[1, 2] <- x[2]
  k0[1, 6] <- x[3]
  k1[1, 6] <- x[4]
  
  # state 2
  k0[2, 3] <- x[5]
  k1[2, 3] <- x[6]
  k0[2, 7] <- x[7]
  k1[2, 7] <- x[8]
  # state 3
  k0[3, 4] <- x[9]
  k1[3, 4] <- x[10]
  # state 4 (from separate estimation)
  k0[4, 5] <- r45[1]
  k1[4, 5] <- r45[2]
  # state 5
  k0[5, 8] <- 0.3
  k1[5, 8] <- 0
  # state 6
  k0[6, 3] <- x[11]
  k1[6, 3] <- x[12]
  k0[6, 4] <- x[13]
  k1[6, 4] <- x[14]
  k0[6, 7] <- x[15]
  k1[6, 7] <- x[16]
  # state 7
  k0[7, 4] <- x[17]
  k1[7, 4] <- x[18]
  
  return(list(k0, k1))
}


TP.f.ATN <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  rec_states <- vector("list", length = (ad.state - 1))
  
  for(i in 1:length(rec_states)){
    # give indices for which states have nonzero transition probabilities
    rec_states[[i]] <- which(k0[i,] > 0)
  }
  
  if(y < int.year) {
    alpha = matrix(1, nrow = (ad.state + 1), ncol = (ad.state + 2)) ##### Ask about this one
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  # we assume rates don't change after age 95
  if(a > 95){
    a = 95
  } else {
    a = a
  }
  p <- matrix(0, (ad.state + 1), (ad.state + 2))
  p[10, 11] <- alpha[10, 11] * 0.167 * (1 - dr.f(a, y, f))   # P(alz -> advanced alz) = 0.167 if you don't die
  p[10, death] <- alpha[10, death] * dr.f(a, y, f)                # P(death) 
  p[11, death] <- alpha[11, death] * (dr.f(a, y, f) + 0.078)     # P(death | advanced alz) = p(death) + 0.078
  for(i in 1:(ad.state - 1)) {
    for(j in rec_states[[i]]){
      # transitions to states other than death
      # age < 65 is default
      p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.f(a, y, f))
      
      # transition to death
      p[i, death] <- alpha[i, death] * dr.f(a, y, f)
      
      # transitions are slightly different if you have mcid
      if(i %in% c(5)){
        p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.f(a, y, f) * mcid))
        p[i, death] <- alpha[i, death] * (dr.f(a, y, f) * mcid)
      }
    }
  }
  for(i in 1:(ad.state + 1)){
    # prob of staying in the same state
    p[i, i] <- 1 - sum(p[i,])
    
    # if(a < 100){
    #   # prob of staying in the same state
    #   p[i, i] <- 1 - sum(p[i,]) }
    # else if (a > 99){ 
    #   #if you're 100 there's no way you're normal
    #   p[1, 1] <- 0
    #   p[i, i] <- 1 - sum(p[i,]) }
  }
  # if you're dead you stay dead
  prob <- rbind(p, c(rep(0, (ad.state + 1)), 1))
  return(prob)
}


##### Now transition probabilities for males
TP.m.ATN <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  rec_states <- vector("list", length = (ad.state - 1))
  
  for(i in 1:length(rec_states)){
    # give indices for which states have nonzero transition probabilities
    rec_states[[i]] <- which(k0[i,] > 0)
  }
  
  if(y < int.year) {
    alpha = matrix(1, nrow = (ad.state + 1), ncol = (ad.state + 2)) ##### Ask about this one
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  if(a > 95){
    a = 95
  } else {
    a = a
  }
  p <- matrix(0, (ad.state + 1), (ad.state + 2))
  p[10, 11] <- alpha[10, 11] * 0.167 * (1 - dr.m(a, y, f))   # P(alz -> advanced alz) = 0.167 if you don't die
  p[10, death] <- alpha[10, death] * dr.m(a, y, f)                # P(death) 
  p[11, death] <- alpha[11, death] * (dr.m(a, y, f) + 0.078)     # P(death | advanced alz) = p(death) + 0.078
  for(i in 1:(ad.state - 1)) {
    for(j in rec_states[[i]]){
      # transitions to states other than death
      p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.m(a, y, f))
      
      # transition to death
      p[i, death] <- alpha[i, death] * dr.m(a, y, f)
      
      # there's a multiplier for P(death) with mcid
      if(i %in% c(5)){
        p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.m(a, y, f) * mcid))
        p[i, death] <- alpha[i, death] * (dr.m(a, y, f) * mcid)
      }
    }
  }
  for(i in 1:(ad.state + 1)){
    # prob of staying in same state
    p[i, i] <- 1 - sum(p[i,])
    
    # if(a < 100){
    #   # prob of staying in the same state
    #   p[i, i] <- 1 - sum(p[i,]) }
    # else if (a > 99){ 
    #   #if you're 100 there's no way you're normal
    #   p[1, 1] <- 0
    #   p[i, i] <- 1 - sum(p[i,]) }
  }
  # if you're dead you stay dead
  prob <- rbind(p, c(rep(0, (ad.state + 1)), 1))
  return(prob)
}

Phi.f.ATN <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  n = a - 30  
  y2 <- y - n - 1
  # prod object will contain the matrix products
  prod <- list()
  prod[[1]] <- TP.f.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
  }
  
  # take first row because they started in state 1 (normal) when they were 30
  phi <- prod[[n + 1]][1, 1:(ad.state + 1)]
  return(phi)
}

# Equal function for males, uses TP.m instead of TP.f
Phi.m.ATN <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  n = a - 30  
  y2 <- y - n - 1
  
  prod <- list()
  prod[[1]] <- TP.m.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
  }
  
  phi <- prod[[n + 1]][1, 1:(ad.state + 1)]
  return(phi)
}



#### These prevalence rates will ONLY be for comparing to the prevalence rates given by Jack (2017)
#### Jack only gives prevalences for the eight states: 1-4 and 6-9
Prevrate.f.ATN <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  phi = Phi.f.ATN(age, y, alpha, int.year, mcid, f, k0, k1)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

Prevrate.m.ATN <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  phi = Phi.m.ATN(age, y, alpha, int.year, mcid, f, k0, k1)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

### These are prevalence functions that will return a matrix if you give it a vector for age
### matrix will be m x n where m = length(age) and n = 8 (for preclinical states)

Prevrate.f.multi.ATN <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, c(1:4,6:9)]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

Prevrate.m.multi.ATN <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, c(1:4,6:9)]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

Prevrate.f.multi.ATN.uncond <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
  }
  
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:9)]
  return(prevrates)
}

Prevrate.m.multi.ATN.uncond <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
  }
  
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:9)]
  return(prevrates)
}

incidence.females <- NULL

incidence.f.multi <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){

  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN(30, y2, intalpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN(30 + i, y2 + i, intalpha, int.year = 2014, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 10] /
                                          sum(prev_curr[1:9])
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  
  return(inc.multi)
  
}

incidence.m.multi <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN(30, y2, intalpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN(30 + i, y2 + i, intalpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 10] /
        sum(prev_curr[1:9])
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  
  return(inc.multi)
  
}



incidence.f <- function(a, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  incidence.females <- (Prevrate.f.ATN(a - 1, y,  alpha, int.year, mcid, f, k0, k1)[5] * 
                          TP.f.ATN(a - 1, y - 1, alpha, int.year, mcid, f, k0, k1)[5, 10]) /
    sum(Prevrate.f.ATN(a - 1, y, alpha, int.year, mcid, f, k0, k1)[1:9])
  return(incidence.females)
}

incidence.males <- NULL
incidence.m <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  incidence.males <- (Prevrate.m.ATN(a - 1, y,  alpha, int.year, mcid, f, k0, k1)[5] * TP.m(a - 1, y - 1, alpha, int.year, mcid, f, k0, k1)[5, 10]) /
    sum(Prevrate.m.ATN(a - 1, y, alpha, int.year, mcid, f, k0, k1)[1:9])
  return(incidence.males)
}

incidence.prevrate.f.uncond <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN(30, y2, intalpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN(30 + i, y2 + i, intalpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 10] /
        sum(prev_curr[1:9])
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:9)]
  return(cbind(prevrates, inc.multi))
  
}

incidence.prevrate.m.uncond <- function(age, y = 2014, alpha = intalpha, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN(30, y2, intalpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN(30 + i, y2 + i, intalpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.m.ATN(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 10] /
        sum(prev_curr[1:9])
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:9)]
  return(cbind(prevrates, inc.multi))
  
}

######### Lifetime!

AR.f.ATN <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  prod <- sum <- vector("list")
  prob.adj.f <- NULL
  pr.rec.adj.f <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    
    prob.adj.f <- rbind(TP.f.ATN(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.f.ATN(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.f)
  }
  
  prod[[1]] <- pr.rec.adj.f(a, y, int.year, alpha, mcid, f, k0, k1)
  for(i in 1:(t - 1)){
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.f(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
  }
  
  sum[[1]] <- prod[[1]]
  for(k in 1:(t - 1)){
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}



# MALES

AR.m.ATN <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  prod <- sum <- vector("list")
  prob.adj <- prob.adj.m <- NULL
  pr.rec.adj.m <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    prob.adj.m <- rbind(TP.m.ATN(a, y, alpha, int.year, mcid, f, k0, k1)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.m.ATN(a, y, alpha, int.year,  mcid, f, k0, k1)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.m)
  }
  for(i in 1:(t - 1)){
    prod[[1]] <- pr.rec.adj.m(a, y, int.year, alpha, mcid, f, k0, k1)
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.m(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
  }
  for(k in 1:(t - 1)){
    sum[[1]] <- prod[[1]]
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}


# The lifetime risk for a person who is age a and is in disease state i at calendar time y, 
# is denoted by LR(a,y,i), and is given by LR(a,y,i)=∑M−an=1ϕ(n|a,y,i). Since the function 
# depends on AR.m and AR.f, it also depends on the effects of interventions (the proportionality 
# constants (or the relative risks)), intervention year int.year, death rate multiplicative 
# factor f and mci.d. The function can also take arguments that specify the transition rates 
# constatnts k0,k1,k012,k112 and k12. This function gives the five year, ten year and lifetime risks.

# vector containing coefficient inputs

# coef.inputs <- c(a12, b12, a16, b16, a18, b18, a23, b23, a34, b34, a45, b45, a510, b510,
#                  a63, b63, a67, b67, a74, b74, a87, b87, a89, b89, a94, b94)

lifetime <- function(age, g, state, int.year = 2014, alpha = intalpha, mcid = 1.65, f = 1, k0, k1){
  
  
  ### Update these in a way that's appropriate
  if(state == "Normal[1]"){
    state = 1
  }else if(state == "Amyloidosis[2]"){
    state = 2
  }else if(state == "A&T[3]"){
    state = 3
  }else if(state == "A&&T&N[4]"){
    state = 4
  }else if(state == "MCI&A&T&N[5]"){
    state = 5
  }else if(state == "Tauopathy[6]"){
    state = 6
  }else if(state == "T&N[7]"){
    state = 7
  }else if(state == "Neurodegeneration[8]"){
    state = 8
  }else if(state == "A&N[9]"){
    state = 9
  }
  
  if(g == "Male"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.m.ATN(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.m.ATN(age, 10, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1,ncol = (ad.state - 1))
    fiveyearrisk <- AR.m.ATN(age, 5, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }else if(g=="Female"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.f.ATN(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.f.ATN(age, 10, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    fiveyearrisk <- AR.f.ATN(age, 5, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }
  return(c(paste(round(lifetimerisk[state], 2), "%"), paste(round(tenyearrisk[state], 2), "%"), 
           paste(round(fiveyearrisk[state], 2), "%")))
}




#################################################

########## REDUCED ATN MODEL 

#################################################

TP.f.ATN.reduced <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  ad.state <- 8
  death <- 10
  rec_states <- vector("list", length = (ad.state - 1))
  
  for(i in 1:length(rec_states)){
    # give indices for which states have nonzero transition probabilities
    rec_states[[i]] <- which(k0[i,] > 0)
  }
  
  if(y < int.year) {
    alpha = matrix(1, nrow = (ad.state + 1), ncol = (ad.state + 2)) ##### Ask about this one
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  # we assume rates don't change after age 95
  if(a > 95){
    a = 95
  } else {
    a = a
  }
  p <- matrix(0, (ad.state + 1), (ad.state + 2))
  p[ad.state, (ad.state + 1)] <- alpha[ad.state, (ad.state + 1)] * 0.167 * (1 - dr.f(a, y, f))   # P(alz -> advanced alz) = 0.167 if you don't die
  p[ad.state, death] <- alpha[ad.state, death] * dr.f(a, y, f)                # P(death) 
  p[(ad.state + 1), death] <- alpha[(ad.state + 1), death] * (dr.f(a, y, f) + 0.078)     # P(death | advanced alz) = p(death) + 0.078
  for(i in 1:(ad.state - 1)) {
    for(j in rec_states[[i]]){
      # transitions to states other than death
      # age < 65 is default
      p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.f(a, y, f))
      
      # transition to death
      p[i, death] <- alpha[i, death] * dr.f(a, y, f)
      
      # transitions are slightly different if you have mcid
      if(i %in% c(5)){
        p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.f(a, y, f) * mcid))
        p[i, death] <- alpha[i, death] * (dr.f(a, y, f) * mcid)
      }
    }
  }
  for(i in 1:(ad.state + 1)){
    # prob of staying in the same state
    p[i, i] <- 1 - sum(p[i,])
    
    # if(a < 100){
    #   # prob of staying in the same state
    #   p[i, i] <- 1 - sum(p[i,]) }
    # else if (a > 99){ 
    #   #if you're 100 there's no way you're normal
    #   p[1, 1] <- 0
    #   p[i, i] <- 1 - sum(p[i,]) }
  }
  # if you're dead you stay dead
  prob <- rbind(p, c(rep(0, (ad.state + 1)), 1))
  return(prob)
}

TP.m.ATN.reduced <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  ad.state <- 8
  death <- 10
  rec_states <- vector("list", length = (ad.state - 1))
  
  for(i in 1:length(rec_states)){
    # give indices for which states have nonzero transition probabilities
    rec_states[[i]] <- which(k0[i,] > 0)
  }
  
  if(y < int.year) {
    alpha = matrix(1, nrow = (ad.state + 1), ncol = (ad.state + 2)) ##### Ask about this one
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  # we assume rates don't change after age 95
  if(a > 95){
    a = 95
  } else {
    a = a
  }
  p <- matrix(0, (ad.state + 1), (ad.state + 2))
  p[ad.state, (ad.state + 1)] <- alpha[ad.state, (ad.state + 1)] * 0.167 * (1 - dr.m(a, y, f))   # P(alz -> advanced alz) = 0.167 if you don't die
  p[ad.state, death] <- alpha[ad.state, death] * dr.m(a, y, f)                # P(death) 
  p[(ad.state + 1), death] <- alpha[(ad.state + 1), death] * (dr.m(a, y, f) + 0.078)     # P(death | advanced alz) = p(death) + 0.078
  for(i in 1:(ad.state - 1)) {
    for(j in rec_states[[i]]){
      # transitions to states other than death
      # age < 65 is default
      p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.m(a, y, f))
      
      # transition to death
      p[i, death] <- alpha[i, death] * dr.m(a, y, f)
      
      # transitions are slightly different if you have mcid
      if(i %in% c(5)){
        p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.m(a, y, f) * mcid))
        p[i, death] <- alpha[i, death] * (dr.m(a, y, f) * mcid)
      }
    }
  }
  for(i in 1:(ad.state + 1)){
    # prob of staying in the same state
    p[i, i] <- 1 - sum(p[i,])
    
    # if(a < 100){
    #   # prob of staying in the same state
    #   p[i, i] <- 1 - sum(p[i,]) }
    # else if (a > 99){ 
    #   #if you're 100 there's no way you're normal
    #   p[1, 1] <- 0
    #   p[i, i] <- 1 - sum(p[i,]) }
  }
  # if you're dead you stay dead
  prob <- rbind(p, c(rep(0, (ad.state + 1)), 1))
  return(prob)
}

Phi.f.ATN.reduced <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  ad.state <- 8
  death <- 10
  n = a - 30  
  y2 <- y - n - 1
  # prod object will contain the matrix products
  prod <- list()
  prod[[1]] <- TP.f.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
  }
  
  # take first row because they started in state 1 (normal) when they were 30
  phi <- prod[[n + 1]][1, 1:(ad.state + 1)]
  return(phi)
}

# Equal function for males, uses TP.m instead of TP.f
Phi.m.ATN.reduced <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  ad.state <- 8
  death <- 10
  n = a - 30  
  y2 <- y - n - 1
  
  prod <- list()
  prod[[1]] <- TP.m.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
  }
  
  phi <- prod[[n + 1]][1, 1:(ad.state + 1)]
  return(phi)
}



#### These prevalence rates will ONLY be for comparing to the prevalence rates given by Jack (2017)
#### Jack only gives prevalences for the eight states: 1-4 and 6-9
Prevrate.f.ATN.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  phi = Phi.f.ATN.reduced(age, y, alpha, int.year, mcid, f, k0, k1)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

Prevrate.m.ATN.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  phi = Phi.m.ATN.reduced(age, y, alpha, int.year, mcid, f, k0, k1)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

### These are prevalence functions that will return a matrix if you give it a vector for age
### matrix will be m x n where m = length(age) and n = 8 (for preclinical states)

Prevrate.f.multi.ATN.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, c(1:4,6:7)]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

Prevrate.m.multi.ATN.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, c(1:4,6:7)]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

Prevrate.f.multi.ATN.uncond.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  ad.state <- 8
  death <- 10
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
  }
  
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:7)]
  return(prevrates)
}

Prevrate.m.multi.ATN.uncond.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  # we go back to when they were 30 and multiply matrices from there
  ad.state <- 8
  death <- 10
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
  }
  
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:7)]
  return(prevrates)
}

incidence.females <- NULL

incidence.f.multi.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  ad.state <- 8
  death <- 10
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN.reduced(30 + i, y2 + i, alpha, int.year = 2014, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN.reduced(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 8] /
        sum(prev_curr[1:7])
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  
  return(inc.multi)
  
}

incidence.m.multi.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  ad.state <- 8
  death <- 10
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN.reduced(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 8] /
        sum(prev_curr[1:7])
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  
  return(inc.multi)
  
}



incidence.f.reduced <- function(a, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  incidence.females <- (Prevrate.f.ATN.reduced(a - 1, y,  alpha, int.year, mcid, f, k0, k1)[5] * 
                          TP.f.ATN.reduced(a - 1, y - 1, alpha, int.year, mcid, f, k0, k1)[5, 8]) /
    sum(Prevrate.f.ATN.reduced(a - 1, y, alpha, int.year, mcid, f, k0, k1)[1:7])
  return(incidence.females)
}

incidence.males <- NULL
incidence.m.reduced <- function(a, y, alpha, int.year, mcid, f, k0, k1){
  incidence.males <- (Prevrate.m.ATN.reduced(a - 1, y,  alpha, int.year, mcid, f, k0, k1)[5] * TP.m.reduced(a - 1, y - 1, alpha, int.year, mcid, f, k0, k1)[5, 8]) /
    sum(Prevrate.m.ATN.reduced(a - 1, y, alpha, int.year, mcid, f, k0, k1)[1:7])
  return(incidence.males)
}

incidence.prevrate.f.uncond.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  ad.state <- 8
  death <- 10
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.f.ATN.reduced(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 8] /
        sum(prev_curr[1:7])
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:7)]
  return(cbind(prevrates, inc.multi))
  
}

incidence.prevrate.m.uncond.reduced <- function(age, y = 2014, alpha = intalpha.reduced, int.year = 2014, mcid = 1.65, f = 1, k0, k1){
  ad.state <- 8
  death <- 10
  # we go back to when they were 30 and multiply matrices from there
  inc.multi <- vector(length = length(age))
  phi.multi <- matrix(nrow = length(age), ncol = (death - 1))
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.ATN.reduced(30, y2, alpha, int.year, mcid, f, k0, k1)
  prev_curr <- vector(length = (death - 1))
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.ATN.reduced(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1)
    if((i + 30 - first) >= 0){
      inc.multi[(i + 30 - first + 1)] <- prev_curr[5] * TP.m.ATN.reduced(30 + i - 1, y2 + i - 1, alpha, int.year, mcid, f, k0, k1)[5, 8] /
        sum(prev_curr[1:7])
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:(death - 1)]
    }
    prev_curr <- prod[[i + 1]][1, 1:(death - 1)] / sum(prod[[i + 1]][1, 1:(death - 1)])
  }
  prevrates <- (phi.multi / rowSums(phi.multi))[, c(1:4, 6:7)]
  return(cbind(prevrates, inc.multi))
  
}

######### Lifetime!

AR.f.ATN.reduced <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  ad.state <- 8
  death <- 10
  prod <- sum <- vector("list")
  prob.adj.f <- NULL
  pr.rec.adj.f <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    
    prob.adj.f <- rbind(TP.f.ATN.reduced(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.f.ATN.reduced(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.f)
  }
  
  prod[[1]] <- pr.rec.adj.f(a, y, int.year, alpha, mcid, f, k0, k1)
  for(i in 1:(t - 1)){
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.f(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
  }
  
  sum[[1]] <- prod[[1]]
  for(k in 1:(t - 1)){
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}



# MALES

AR.m.ATN.reduced <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  prod <- sum <- vector("list")
  prob.adj <- prob.adj.m <- NULL
  pr.rec.adj.m <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    prob.adj.m <- rbind(TP.m.ATN.reduced(a, y, alpha, int.year, mcid, f, k0, k1)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.m.ATN.reduced(a, y, alpha, int.year,  mcid, f, k0, k1)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.m)
  }
  for(i in 1:(t - 1)){
    prod[[1]] <- pr.rec.adj.m(a, y, int.year, alpha, mcid, f, k0, k1)
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.m(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
  }
  for(k in 1:(t - 1)){
    sum[[1]] <- prod[[1]]
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}


# The lifetime risk for a person who is age a and is in disease state i at calendar time y, 
# is denoted by LR(a,y,i), and is given by LR(a,y,i)=∑M−an=1ϕ(n|a,y,i). Since the function 
# depends on AR.m and AR.f, it also depends on the effects of interventions (the proportionality 
# constants (or the relative risks)), intervention year int.year, death rate multiplicative 
# factor f and mci.d. The function can also take arguments that specify the transition rates 
# constatnts k0,k1,k012,k112 and k12. This function gives the five year, ten year and lifetime risks.

# vector containing coefficient inputs

# coef.inputs <- c(a12, b12, a16, b16, a18, b18, a23, b23, a34, b34, a45, b45, a510, b510,
#                  a63, b63, a67, b67, a74, b74, a87, b87, a89, b89, a94, b94)

lifetime.reduced <- function(age, g, state, int.year = 2014, alpha = intalpha.reduced, mcid = 1.65, f = 1, k0, k1){
  
  
  ### Update these in a way that's appropriate
  if(state == "Normal[1]"){
    state = 1
  }else if(state == "Amyloidosis[2]"){
    state = 2
  }else if(state == "A&T[3]"){
    state = 3
  }else if(state == "A&&T&N[4]"){
    state = 4
  }else if(state == "MCI&A&T&N[5]"){
    state = 5
  }else if(state == "Tauopathy[6]"){
    state = 6
  }else if(state == "T&N[7]"){
    state = 7
  }else if(state == "Neurodegeneration[8]"){
    state = 8
  }else if(state == "A&N[9]"){
    state = 9
  }
  
  if(g == "Male"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.m.ATN.reduced(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.m.ATN.reduced(age, 10, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1,ncol = (ad.state - 1))
    fiveyearrisk <- AR.m.ATN.reduced(age, 5, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }else if(g=="Female"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.f.ATN.reduced(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.f.ATN.reduced(age, 10, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    fiveyearrisk <- AR.f.ATN.reduced(age, 5, 2017, int.year, alpha, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }
  return(c(paste(round(lifetimerisk[state], 2), "%"), paste(round(tenyearrisk[state], 2), "%"), 
           paste(round(fiveyearrisk[state], 2), "%")))
}







########################################

####### AN MODEL with Nada's old code

########################################

TP.f.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  
  if(y < int.year){
    alpha = matrix(1, nrow = 8, ncol = 9)
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  if(a > 109){
    a = 110
  } else {
    a = a
  }
  p <- matrix(0, 8, 9)
  p[7, 8] <- alpha[7, 8] * 0.167 * (1 - dr.f(a, y, f))
  p[7, 9] <- alpha[7, 9] * dr.f(a, y, f)
  p[8, 9] <- alpha[8, 9] * (dr.f(a, y, f) + 0.078)
  # p[7,8] <- 0
  # p[7,9] <- 1
  
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.f(a, y, f))}
        else if(65 <= a & a <= 75){
          p[1, 2] <- alpha[1, 2] * k012 * exp(k112 * a) * (1 - dr.f(a, y, f))
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.f(a, y, f))
        }
        else if(a > 75){
          p[1, 2] <- alpha[1, 2] * k12 * (1 - dr.f(a, y, f))
          p[3, 6] <- alpha[3, 6] * k0[3, 6] * exp(k1[3, 6] * 75) * (1 - dr.f(a, y, f))
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.f(a, y, f))
        }
        p[i, 9] <- alpha[i, 9] * dr.f(a, y, f)
        if(i %in% 5:6){
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.f(a, y, f) * mcid))
          p[i, 9] <- alpha[i, 9] * (dr.f(a, y, f) * mcid)
        }
      }
    }
  }
  for(i in 1:8){
    if(a < 100){
      p[i, i] <- 1 - sum(p[i, 1:9]) }
    else if (a > 99){
      p[1, 1] <- 0
      p[i, i] <- 1 - sum(p[i, 1:9]) }
  }
  prob <- rbind(p, c(rep(0, 8), 1))
  return(prob)
}
TP.m.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  if(y < int.year){
    alpha = matrix(1, nrow = 8, ncol = 9)
  } else if(y >= int.year){
    alpha = alpha
  }
  if(y > 2014){
    y = 2014
  } else if (y < 1933){
    y <- 1933
  }
  else {y = y}
  if(a > 109){
    a = 110
  } else {
    a = a
  }
  p <- matrix(0, 8, 9)
  p[7, 8] <- alpha[7, 8] * 0.167 * (1 - dr.m(a, y, f))
  p[7, 9] <- alpha[7, 9] * dr.m(a, y, f)
  p[8, 9] <- alpha[8, 9] * (dr.m(a, y, f) + 0.078)
  # p[7, 9] <- 1
  # p[7, 8] <- 0
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.m(a, y, f))}
        else if(65 <= a & a <= 75){
          p[1, 2] <- alpha[1, 2] * k012 * exp(k112 * a) * (1 - dr.m(a, y, f))
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.m(a, y, f))
        }
        else if(a>75){
          p[1, 2] <- alpha[1, 2] * 0.07 * (1 - dr.m(a, y, f))
          p[3, 6] <- alpha[3, 6] * k0[3, 6] * exp(k1[3, 6] * 75) * (1 - dr.m(a, y, f))
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - dr.m(a, y, f))
        }
        p[i, 9] <- alpha[i, 9] * dr.m(a, y, f)
        if(i %in% 5:6){
          p[i, j] <- alpha[i, j] * k0[i, j] * exp(k1[i, j] * a) * (1 - (dr.m(a, y, f) * mcid))
          p[i, 9] <- alpha[i, 9] * (dr.m(a, y, f) * mcid)
        }
      }
    }
  }
  for(i in 1:8){
    if(a < 100){
      p[i, i] <- 1 - sum(p[i, 1:9]) }
    else if (a > 99){
      p[1, 1] <- 0
      p[i, i] <- 1 - sum(p[i, 1:9]) }
  }
  prob <- rbind(p, c(rep(0, 8), 1))
  return(prob)
}

#### These aren't really used but here they are
Phi.f.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  prod <- list()
  n = a - 30  
  y2 <- y - n - 1
  for(i in 1:(n)){
    prod[[1]] <- TP.f.AN(30, y2, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.AN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)}
  phi <- prod[[n + 1]][1, 1:8]
  return(phi)
}

Phi.m.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  prod <- list()
  n = a - 30  
  y2 <- y - n - 1
  for(i in 1:(n)){
    prod[[1]] <- TP.m.AN(30, y2, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.AN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)}
  phi <- prod[[n + 1]][1, 1:8]
  return(phi)
}

Prevrate.f.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  phi = Phi.f.AN(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

Prevrate.m.AN <- function(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  phi = Phi.m.AN(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
  prevalence <- phi / sum(phi)
  return(c(prevalence))
}

### Prevalences for all states up through advanced AD dementia
### give an age vector and it'll give prevalences for states 1-8 for each age

Prevrate.f.multi.AN <- function(age, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  prod <- list()
  # prod object will contain the matrix products
  prod[[1]] <- TP.f.AN(30, y2, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.f.AN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:8]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

Prevrate.m.multi.AN <- function(age, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12){
  # we go back to when they were 30 and multiply matrices from there
  phi.multi <- matrix(nrow = length(age), ncol = 8)
  n = max(age) - 30  
  first = min(age)
  y2 <- y - n - 1
  # prod object will contain the matrix products
  prod[[1]] <- TP.m.AN(30, y2, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
  for(i in 1:(n)){
    prod[[i + 1]] <- (prod[[i]]) %*% TP.m.AN(30 + i, y2 + i, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)
    if((i + 30 - first) >= 0)
      phi.multi[(i + 30 - first + 1),] <- prod[[i + 1]][1, 1:8]
  }
  
  prevrates <- phi.multi / rowSums(phi.multi)
  return(prevrates)
}

AR.f.AN <- function(a, t, y, int.year, alpha, mcid = 1.65, f = 1, k0, k1, k012, k112, k12){
  ad.state <- 7
  
  prod <- sum <- vector("list")
  prob.adj.f <- NULL
  pr.rec.adj.f <- function(a, y, int.year, alpha, mcid, f, k0, k1, k012, k112, k12){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    
    prob.adj.f <- rbind(TP.f.AN(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.f.AN(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.f)
  }
  
  prod[[1]] <- pr.rec.adj.f(a, y, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)
  for(i in 1:(t - 1)){
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.f(a + i, y + i, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)
  }
  
  sum[[1]] <- prod[[1]]
  for(k in 1:(t - 1)){
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}



# MALES

AR.m.AN <- function(a, t, y, int.year, alpha, mcid, f, k0, k1, k012, k112, k12){
  ad.state <- 7
  prod <- sum <- vector("list")
  prob.adj <- prob.adj.m <- NULL
  pr.rec.adj.m <- function(a, y, int.year, alpha, mcid, f, k0, k1, k012, k112, k12){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    prob.adj.m <- rbind(TP.m.AN(a, y, alpha, int.year, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.m.AN(a, y, alpha, int.year,  mcid, f, k0, k1, k012, k112, k12)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj.m)
  }
  for(i in 1:(t - 1)){
    prod[[1]] <- pr.rec.adj.m(a, y, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj.m(a + i, y + i, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)
  }
  for(k in 1:(t - 1)){
    sum[[1]] <- prod[[1]]
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}


# The lifetime risk for a person who is age a and is in disease state i at calendar time y, 
# is denoted by LR(a,y,i), and is given by LR(a,y,i)=∑M−an=1ϕ(n|a,y,i). Since the function 
# depends on AR.m and AR.f, it also depends on the effects of interventions (the proportionality 
# constants (or the relative risks)), intervention year int.year, death rate multiplicative 
# factor f and mci.d. The function can also take arguments that specify the transition rates 
# constatnts k0,k1,k012,k112 and k12. This function gives the five year, ten year and lifetime risks.

# vector containing coefficient inputs

# coef.inputs <- c(a12, b12, a16, b16, a18, b18, a23, b23, a34, b34, a45, b45, a510, b510,
#                  a63, b63, a67, b67, a74, b74, a87, b87, a89, b89, a94, b94)

lifetime.AN <- function(age, g, state, int.year, alpha = intalpha, mcid = 1.65, f = 1, 
                        k0, k1, k012, k112, k12){
  
  ad.state <- 7
  ### Update these in a way that's appropriate
  if(state == "Normal[1]"){
    state = 1
  }else if(state == "Amyloidosis[2]"){
    state = 2
  }else if(state == "Neurodegeneration[3]"){
    state = 3
  }else if(state == "A&N[4]"){
    state = 4
  }else if(state == "MCI&A&N[5]"){
    state = 5
  }else if(state == "MCI&N[6]"){
    state = 6
  }
  
  if(g == "Male"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.m.AN(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.m.AN(age, 10, 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1,ncol = (ad.state - 1))
    fiveyearrisk <- AR.m.AN(age, 5, 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
  }else if(g == "Female"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.f.AN(age, (110 - age), 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.f.AN(age, 10, 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    fiveyearrisk <- AR.f.AN(age, 5, 2017, int.year, alpha, mcid, f, k0, k1, k012, k112, k12)[1:(ad.state - 1), ad.state] * 100
  }
  return(c(paste(round(lifetimerisk[state], 2), "%"), paste(round(tenyearrisk[state], 2), "%"), 
           paste(round(fiveyearrisk[state], 2), "%")))
}

k0.AN<-matrix(0,6,7)
k1.AN<-matrix(0,6,7)
k0.AN[1,2]<-0.000149
k1.AN[1,2]<-0.086
k0.AN[1,3]<-7.531246e-06
k1.AN[1,3]<-0.117866
k0.AN[3,4]<-0.00012646
k1.AN[3,4]<-0.081938
k0.AN[2,4]<-0.001109085
k1.AN[2,4]<-0.0616514
k0.AN[3,6]<-0.003276065
k1.AN[3,6]<-0.01981314
k0.AN[6,7]<-0.09
k1.AN[6,7]<-0
k0.AN[4,5]<-2.595831e-05
k1.AN[4,5]<-0.096183077
k0.AN[5,7]<-0.3
k1.AN[5,7]<-0
k012.AN<-0.00105
k112.AN<-0.05596
k12.AN<-0.07


