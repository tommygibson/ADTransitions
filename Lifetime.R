######## LIFETIME RISKS IN ATN MODEL

source("BrookFuncs.R")

readRDS("weight.1.opt_04.21.2021.rds")

AR.f.ATN <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  prod <- sum <- vector("list")
  prob.adj <- prob.adj.m <- NULL
  pr.rec.adj <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    
    prob.adj <- rbind(TP.f.ATN(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[1:(ad.state - 1),], 
                      c(rep(0, (ad.state + 1)), 1), 
                      TP.f.ATN(a, y, alpha = alpha, int.year, mcid = 1.65, f = 1, k0, k1)[(ad.state + 1),], 
                      c(rep(0, (ad.state + 1)), 1))
    return(prob.adj)
  }
  
  prod[[1]] <- pr.rec.adj(a, y, int.year, alpha, mcid, f, k0, k1)
  for(i in 1:(t - 1)){
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
  }
  
  sum[[1]] <- prod[[1]]
  for(k in 1:(t - 1)){
    sum[[k + 1]] <- sum[[k]] + prod[[k + 1]]
  }
  
  return(sum[[t]])
}



# MALES

AR.m.ATN <- function(a, t, y, int.year, alpha, mcid, f, k0, k1){
  prob.adj <- prob.adj.m <- NULL
  pr.rec.adj.males <- function(a, y, int.year, alpha, mcid, f, k0, k1){
    
    # here we set the 10th row (transition from AD dementia) to immediately transition them to death
    prob.adj.m <- rbind(TP.m.ATN(a, y, alpha, int.year, mcid, f, k0, k1)[1:(ad.state - 1),], 
                        c(rep(0, (ad.state + 1)), 1), 
                        TP.m.ATN(a, y, alpha, int.year,  mcid, f, k0, k1)[(ad.state + 1),], 
                        c(rep(0, (ad.state + 1)), 1))
    return(prob.adj)
  }
  for(i in 1:(t - 1)){
    prod[[1]] <- pr.rec.adj(a, y, int.year, alpha, mcid, f, k0, k1)
    prod[[i + 1]] <- (prod[[i]]) %*% pr.rec.adj(a + i, y + i, int.year, alpha, mcid, f, k0, k1)
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

lifetime <- function(age, g, state, int.year, alpha.inputs, mcid, f, k0, k1){

  
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
    lifetimerisk <- AR.m.ATN(age, (110 - age), 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.m.ATN(age, 10, 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1,ncol = (ad.state - 1))
    fiveyearrisk <- AR.m.ATN(age, 5, 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }else if(g=="Female"){
    lifetimerisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    lifetimerisk <- AR.f.ATN(age, (110 - age), 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    tenyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    tenyearrisk <- AR.f.ATN(age, 10, 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
    fiveyearrisk <- matrix(0, nrow = 1, ncol = (ad.state - 1))
    fiveyearrisk <- AR.f.ATN(age, 5, 2017, int.year, alpha.inputs, mcid, f, k0, k1)[1:(ad.state - 1), ad.state] * 100
  }
  return(c(paste(round(lifetimerisk[state], 2), "%"), paste(round(tenyearrisk[state], 2), "%"), 
           paste(round(fiveyearrisk[state], 2), "%")))
}
