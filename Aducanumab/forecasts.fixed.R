library(ggplot2)
library(mgcv)
library(dplyr)
library(magrittr)

rates<-read.csv("U.S.mortality.rates.csv")
male.census<-read.csv("male.census.csv")
female.census<-read.csv("female.census.csv")

#females
dr.f<-function(a,y,f){
  rates.sub<-subset(rates, rates[,1]==y)
  d.f<-f*rates.sub[a+1,3]
  return(d.f)
}
#males
dr.m<-function(a,y,f){
  rates.sub<-subset(rates, rates[,1]==y)
  d.f<-f*rates.sub[a+1,4]
  return(d.f)
}

#posible interventions
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
intalpha=alpha.nine(1,1,1,1,1,1,1,1,1)

prob<-NULL
#constants of the probabilities of transition from state i to j
k0<-matrix(0,6,7)
k1<-matrix(0,6,7)
k0[1,2]<-0.000149
k1[1,2]<-0.086
k0[1,3]<-7.531246e-06
k1[1,3]<-0.117866
k0[3,4]<-0.00012646
k1[3,4]<-0.081938
k0[2,4]<-0.001109085
k1[2,4]<-0.0616514
k0[3,6]<-0.003276065
k1[3,6]<-0.01981314
k0[6,7]<-0.09
k1[6,7]<-0
k0[4,5]<-2.595831e-05
k1[4,5]<-0.096183077
k0[5,7]<-0.3
k1[5,7]<-0
k012<-0.00105
k112<-0.05596
k12<-0.07
TP.f<-function(a,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  if(y<int.year){
    alpha=matrix(1,nrow=8,ncol=9)
  }else if(y>=int.year){
    alpha=alpha
  }
  if(y>2014){
    y=2014
  } else if (y<1933){
    y<-1933
  }
  else {y=y}
  if(a>109){
    a=110
  } else {
    a=a
  }
  p<-matrix(0,8,9)
  p[7,8]<-alpha[7,8]*0.167*(1-dr.f(a,y,f))
  p[7,9]<-alpha[7,9]*dr.f(a,y,f)
  p[8,9]<-alpha[8,9]*(dr.f(a,y,f)+0.078)
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y,f))}
        else if(65 <= a & a <= 75){
          p[1,2]<-alpha[1,2]*k012*exp(k112*a)*(1-dr.f(a,y,f))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y,f))
        }
        else if(a>75){
          p[1,2]<-alpha[1,2]*k12*(1-dr.f(a,y,f))
          p[3,6]<-alpha[3,6]*k0[3,6]*exp(k1[3,6]*75)*(1-dr.f(a,y,f))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.f(a,y,f))
        }
        p[i,9]<-alpha[i,9]*dr.f(a,y,f)
        if(i %in% 5:6){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-(dr.f(a,y,f)*mcid))
          p[i,9]<-alpha[i,9]*(dr.f(a,y,f)*mcid)
        }
      }
    }
  }
  for(i in 1:8){
    if(a<100){
      p[i,i]<-1-sum(p[i, 1:9]) }
    else if (a>99){
      p[1,1]<-0
      p[i,i]<-1-sum(p[i, 1:9]) }
  }
  prob<-rbind(p,c(rep(0,8),1))
  return(prob)
}
TP.m<-function(a,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  if(y<int.year){
    alpha=matrix(1,nrow=8,ncol=9)
  }else if(y>=int.year){
    alpha=alpha
  }
  if(y>2014){
    y=2014
  } else if (y<1933){
    y<-1933
  }
  else {y=y}
  if(a>109){
    a=110
  } else {
    a=a
  }
  p<-matrix(0,8,9)
  p[7,8]<-alpha[7,8]*0.167*(1-dr.m(a,y,f))
  p[7,9]<-alpha[7,9]*dr.m(a,y,f)
  p[8,9]<-alpha[8,9]*(dr.m(a,y,f)+0.078)
  for(i in 1:6){
    for(j in 1:7){
      if(i < j){
        if(a < 65){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y,f))}
        else if(65 <= a & a <= 75){
          p[1,2]<-alpha[1,2]*k012*exp(k112*a)*(1-dr.m(a,y,f))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y,f))
        }
        else if(a>75){
          p[1,2]<-alpha[1,2]*0.07*(1-dr.m(a,y,f))
          p[3,6]<-alpha[3,6]*k0[3,6]*exp(k1[3,6]*75)*(1-dr.m(a,y,f))
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-dr.m(a,y,f))
        }
        p[i,9]<-alpha[i,9]*dr.m(a,y,f)
        if(i %in% 5:6){
          p[i,j]<-alpha[i,j]*k0[i,j]*exp(k1[i,j]*a)*(1-(dr.m(a,y,f)*mcid))
          p[i,9]<-alpha[i,9]*(dr.m(a,y,f)*mcid)
        }
      }
    }
  }
  for(i in 1:8){
    if(a<100){
      p[i,i]<-1-sum(p[i, 1:9]) }
    else if (a>99){
      p[1,1]<-0
      p[i,i]<-1-sum(p[i, 1:9]) }
  }
  prob<-rbind(p,c(rep(0,8),1))
  return(prob)
}

Phi.f<-function(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  n=a-30  
  y2<-y-n-1
  for(i in 1:(n)){
    prod[[1]]<-TP.f(30,y2,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
    prod[[i+1]]<-(prod[[i]])%*%TP.f(30+i, y2+i,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)}
  phi<-prod[[n+1]][1,1:8]
  return(phi)
}

Prevrate.f<-function(a,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  phi=Phi.f(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
  prevalence<-phi/sum(phi)
  return(c(prevalence))
}
#males
Phi.m<-function(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  n=a-30  
  y2<-y-n-1
  for(i in 1:(n)){
    prod[[1]]<-TP.m(30,y2,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
    prod[[i+1]]<-(prod[[i]])%*%TP.m(30+i, y2+i,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)}
  phi<-prod[[n+1]][1,1:8]
  return(phi)
}
Prevrate.m<-function(a,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  phi=Phi.m(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
  prevalence<-phi/sum(phi)
  return(c(prevalence))
}

#females
incidence.females<-NULL
Incidence.f<-function(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  incidence.females<-((Prevrate.f(a-1,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6]*TP.f(a-1,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6,7])+
                        (Prevrate.f( a-1,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5]*TP.f(a-1,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5,7]))/
    sum(Prevrate.f(a-1,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[1:6])
  return(incidence.females)
}
incidence.males<-NULL
Incidence.m<-function(a,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  
  incidence.males<-((Prevrate.m(a-1,y, alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6]*TP.m(a-1,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6,7])+
                      (Prevrate.m( a-1,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5]*TP.m(a-1,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5,7]))/
    sum(Prevrate.m(a-1,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)[1:6])
  return(incidence.males)
}


#The number of people in the U.S at age a, in year y
Prevnum.f<-function(a, y, alpha,int.year,mcid,f,female.census,k0,k1,k012,k112,k12){
  pf<-num.females<-prevalencef<-vector("list") 
  for(k in 1:(a-30)){
    pf[[1]]<-((TP.f(30,y-(a-30)-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    pf[[k+1]]<-(pf[[k]]%*%(TP.f(30+k,(y-(a-30)+k-1),alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalencef[[1]]<-pf[[1]][1,1:8]/sum(pf[[(1)]][1,1:8])
    prevalencef[[k+1]]<-pf[[k+1]][1,1:8]/sum(pf[[(k+1)]][1,1:8])
    num.females[[1]]<-prevalencef[[1]]*female.census[(y-2014+1),(3+(30))]
    num.females[[k+1]]<-prevalencef[[k]]*female.census[(y-2014+1),(3+(30+k))]}
  return(num.females)
}

Prevnum.m<-function(a, y, alpha,int.year,mcid,f,male.census,k0,k1,k012,k112,k12){
  pm<-num.males<-prevalencem<-vector("list") 
  for(k in 1:(a-30)){
    pm[[1]]<-((TP.m(30,y-(a-30)-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    pm[[k+1]]<-(pm[[k]]%*%(TP.m(30+k,(y-(a-30)+k-1),alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalencem[[1]]<-pm[[1]][1,1:8]/sum(pm[[(1)]][1,1:8])
    prevalencem[[k+1]]<-pm[[k+1]][1,1:8]/sum(pm[[(k+1)]][1,1:8])
    num.males[[1]]<-prevalencem[[1]]*male.census[(y-2014+1),(3+(30))]
    num.males[[k+1]]<-prevalencem[[k]]*male.census[(y-2014+1),(3+(30+k))]}
  return(num.males)
}

Prevnum.Total<-function(prob1,prob2, y, alpha,i,trans6to7.f,trans5to7.f,trans6to7.m,
                        trans5to7.m,int.year,mcid,f,female.census,male.census,k0,k1,k012,k112,k12){
  a=109
  p<-pf<-num.males<-num.females<-num.tot<-prevalencef<- prevalence<-inc.males<-sum.inc.males<-
    mci4<-inc.females<-sum.inc.females<-inc.tot<-inc.mci4<-vector("list")   
  for(k in 1:(a-30)){
    p[[1]]<-((TP.m(30,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalence[[1]]<-p[[1]][1,1:8]/sum(p[[(1)]][1,1:8])
    num.males[[1]]<-prevalence[[1]]*male.census[(i+1),(3+(30))]
    inc.males[[1]]<-0
    p[[k+1]]<-(prob1[[k]]%*%(TP.m(29+k,(y-1),alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalence[[k+1]]<-p[[k+1]][1,1:8]/sum(p[[(k+1)]][1,1:8])
    #The number of persons in the United States by disease state at age a in year y, gender g
    num.males[[k+1]]<-prevalence[[k]]*male.census[(i+1),(3+(30+k))]
    inc.males[[k+1]]<-(num.males[[k+1]][6]*trans6to7.m[i,k+1])+
      (num.males[[k+1]][5]*trans5to7.m[i,k+1])
    pf[[1]]<-((TP.f(30,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalencef[[1]]<-pf[[1]][1,1:8]/sum(pf[[(1)]][1,1:8])
    num.females[[1]]<-prevalencef[[1]]*female.census[(i+1),(3+(30))]
    inc.females[[1]]<-0
    pf[[k+1]]<-(prob2[[k]]%*%(TP.f(29+k,(y-1),alpha,int.year,mcid,f,k0,k1,k012,k112,k12)))
    prevalencef[[k+1]]<-pf[[k+1]][1,1:8]/sum(pf[[(k+1)]][1,1:8])
    num.females[[k+1]]<-prevalencef[[k]]*female.census[(i+1),(3+(30+k))]
    #The incidence of  Alzheimer's disease at age a in year y (numbers of people, in millions)   
    inc.females[[k+1]]<-(num.females[[k+1]][6]*trans6to7.f[i,k+1])+
      (num.females[[k+1]][5]*trans5to7.f[i,k+1])
    mci4[[1]]<-0
    mci4[[k+1]]<-((num.females[[k+1]][5]*trans5to7.f[i,k+1])+
                    (num.males[[k+1]][5]*trans5to7.m[i,k+1]))
  }
  sum.females<-Reduce('+', num.females)
  sum.males<-Reduce('+', num.males) 
  sum.tot<-sum.males+sum.females
  sum.inc.males<-Reduce('+', inc.males)
  sum.inc.females<-Reduce('+', inc.females)
  inc.tot<-sum.inc.males+sum.inc.females
  inc.mci4<-Reduce('+', mci4)
  return(c(p,pf,sum.tot,inc.tot, inc.mci4, num.females,num.males))
}  


TPN.m<-function(a,n, y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  prod<-vector("list")
  for(i in 1:(n)){
    prod[[1]]<-TP.m(a,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
    prod[[i+1]]<-(prod[[i]])%*%TP.m(a+i, y+i-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
  }
  return(prod)
}

TPN.f<-function(a,n,y,alpha,int.year,mcid,f,k0,k1,k012,k112,k12){
  prod<-vector("list")
  for(i in 1:(n)){
    prod[[1]]<-TP.f(a,y-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
    prod[[i+1]]<-(prod[[i]])%*%TP.f(a+i, y+i-1,alpha,int.year,mcid,f,k0,k1,k012,k112,k12)
  }
  return(prod)
}
#define starting year and ending year for projections
year=2014
# 
# year2=2016
# prob1<-prob2<-vector("list")
# for(j in 1:(109-30)){
#   prob1[[1]]<-TP.m(30,year-1,intalpha, 2017, 1.65, 1, k0, k1, k012, k112, k12)
#   prob2[[1]]<-TP.f(30,year-1,intalpha, 2017, 1.65, 1, k0, k1, k012, k112, k12)
#   prob1[[j+1]]<-TPN.m(a=30,n=j,(year-j),alpha=alpha.nine(1,1,1,1,1,1,1,1,1), 2018, 1.65, 1, k0, k1, k012, k112, k12)[[j+1]]
#   prob2[[j+1]]<-TPN.f(a=30,n=j,(year-j),alpha=alpha.nine(1,1,1,1,1,1,1,1,1), 2018, 1.65, 1, k0, k1, k012, k112, k12)[[j+1]]
# }
# for (i in 1:80){
#   write.csv(prob1[[i]],paste("prob1",i,".csv"), row.names = FALSE)
#   write.csv(prob2[[i]],paste("prob2",i,".csv"), row.names = FALSE)
# }
read.csv("prob2 1 .csv")

probsaved<-list.files(path = here("probs", pattern=".csv$"))
split <- strsplit(probsaved, "prob1 ") 

split <- as.numeric(sapply(split, function(x) x <- sub(".csv", "", x[2])))
probsaved1.order <- probsaved[order(split)]
split2 <- strsplit(probsaved, "prob2 ") 

split2 <- as.numeric(sapply(split2, function(x) x <- sub(".csv", "", x[2])))
probsaved2.order <- probsaved[order(split2)]

prob1s<-prob2s<-prob1<-prob2<-vector("list")
for(i in 1:80){
  prob1s[[i]]<-as.matrix(read.csv(file=here("probs", probsaved1.order[i])))
  colnames(prob1s[[i]]) <- NULL
  prob2s[[i]]<-as.matrix(read.csv(file=here("probs", probsaved2.order[i])))
  colnames(prob2s[[i]]) <- NULL
  prob1s[[i]]<-data.matrix(prob1s[[i]])
  prob2s[[i]]<-data.matrix(prob2s[[i]])
}

prob1[[1]]<-prob1s
prob2[[1]]<-prob2s

prevnum.table<-function(int.year, year2,one.two,  two.four, four.five, one.three, three.four,three.six, 
                        five.seven, six.seven, seven.eight,mcid,f,female.census,male.census,a12,a13,a34,a24,
                        a36,a67,a45,a57,b12,b13,b34,b24,
                        b36,b67,b45,b57,k012,k112,k12){
  k0<-matrix(0,6,7)
  k1<-matrix(0,6,7)
  k0[1,2]<-a12
  k1[1,2]<-b12
  k0[1,3]<-a13
  k1[1,3]<-b13
  k0[3,4]<-a34
  k1[3,4]<-b34
  k0[2,4]<-a24
  k1[2,4]<-b24
  k0[3,6]<-a36
  k1[3,6]<-b36
  k0[6,7]<-a67
  k1[6,7]<-b67
  k0[4,5]<-a45
  k1[4,5]<-b45
  k0[5,7]<-a57
  k1[5,7]<-b57
  k012<-0.00105
  k112<-0.05596
  k12<-0.07
  intalpha<-matrix(1,nrow=8,ncol=9)
  intalpha[1,2]<-one.two
  intalpha[2,4]<-two.four
  intalpha[4,5]<-four.five
  intalpha[1,3]<-one.three
  intalpha[3,4]<-three.four
  intalpha[3,6]<-three.six
  intalpha[5,7]<-five.seven
  intalpha[6,7]<-six.seven
  intalpha[7,8]<-seven.eight
  
  trans6to7.m<-matrix(0,year2-year,(109-29))
  trans5to7.m<-matrix(0,year2-year,(109-29))
  trans6to7.f<-matrix(0,year2-year,(109-29))
  trans5to7.f<-matrix(0,year2-year,(109-29))
  fit<-sum.year<-inc.year<-malesnum<-femalesnum<-vector("list")
  for(l in 1:(109-29)){
    trans6to7.f[1:(year2-year),l]<-TP.f(28+l, 2018, intalpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6,7]
    trans5to7.f[1:(year2-year),l]<-TP.f(28+l, 2018, intalpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5,7]
    trans6to7.m[1:(year2-year),l]<-TP.m(28+l, 2018, intalpha,int.year,mcid,f,k0,k1,k012,k112,k12)[6,7]
    trans5to7.m[1:(year2-year),l]<-TP.m(28+l, 2018, intalpha,int.year,mcid,f,k0,k1,k012,k112,k12)[5,7]
  }
  
  for(i in 1:(year2-year)){
    fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha,i=i,
                            trans6to7.f,trans5to7.f,trans6to7.m,trans5to7.m,int.year,mcid,f,female.census,male.census,
                            k0,k1,k012,k112,k12)
    prob1[[i+1]]<-fit[[i]][1:80]
    prob2[[i+1]]<-fit[[i]][81:160]
    sum.year[[i]]<-fit[[i]][161:170]
    femalesnum[[i]]<-fit[[i]][171:250]
    malesnum[[i]]<-fit[[i]][251:330]
  }
  sum<-unlist(sum.year)
  yearvar<-as.character((year+1):year2)
  table<-data.frame(cbind((yearvar),round(matrix((sum)/1000000, nrow=year2-year, ncol=10,byrow=T),2)))
  colnames(table)<-(c("Year","Normal[1]", "A[2]","N[3]","A & N[4]","MCI & A & N[5]",
                      "MCI & N[6]","Early AD[7]","Late AD[8]","Total Incidence",
                      "Incidence from State 5"))
  colnames(table) <- gsub(".", " ", colnames(table), fixed=TRUE)
  #rownames(table)[1:(year2-year)]<-paste("Year",(year+1):year2,sep="")
  tabout<-table[-(1:(int.year-2015)),]
  femalesyear=matrix(0,nrow=80, ncol=8)
  for(i in 1:80){
    femalesyear[i,]<-femalesnum[[(int.year-2017+1)]][[i]]
  }
  malesyear=matrix(0,nrow=80, ncol=8)
  for(i in 1:79){
    malesyear[i,]<-malesnum[[(int.year-2017+1)]][[i]]
  }
  totalyear=malesyear+femalesyear
  plot<-data.frame(cbind(age=c(seq(30,109,1)),totalyear/1000000))
  
  
  totalpooled=cbind(totalyear[,1:4],totalyear[,5]+totalyear[,6],totalyear[,7:8])
  plotpool<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80), rep(7,80)),c(totalpooled[,2:7])/1000000))
  dpool=data.frame(lt=c(rep("solid",5),rep( "dashed", 1)))
  plotyear<-ggplot(plotpool, aes(age,V3+0.01, colour=factor(group),linetype=factor(group))) +
    geom_smooth(method="gam", fill=NA, formula=y~s(x), size=0.92)+
    #scale_y_log10(breaks =c(50,100,200,400,1000),
    #            labels = c(50,100,200,400,1000))+
    xlab("Age")+
    ylab(paste("U.S. Prevalence (millions) in",int.year))+
    scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
    scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
    scale_color_manual(name="",values=c("#33CC33","#0000CC","#3399CC","#000000","#FF66CC","#CC0033"),
                       labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
    scale_linetype_manual(values=c(rep("solid",3),rep( "dotted", 1),rep("solid",2)),name="",
                          labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
    theme(legend.key = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key.size = unit(1.5, 'lines'),
          legend.text=element_text(size=13),text = element_text(size=12),axis.text=element_text(color = "black", size=12))
  
  return(list(tabout, plotyear))
}



prevnum.table.fixed <- function(int.year, year2,one.two,  two.four, four.five, one.three, three.four,three.six, 
                        five.seven, six.seven, seven.eight,mcid,f,female.census,male.census,a12,a13,a34,a24,
                        a36,a67,a45,a57,b12,b13,b34,b24,
                        b36,b67,b45,b57,k012,k112,k12){
  k0<-matrix(0,6,7)
  k1<-matrix(0,6,7)
  k0[1,2]<-a12
  k1[1,2]<-b12
  k0[1,3]<-a13
  k1[1,3]<-b13
  k0[3,4]<-a34
  k1[3,4]<-b34
  k0[2,4]<-a24
  k1[2,4]<-b24
  k0[3,6]<-a36
  k1[3,6]<-b36
  k0[6,7]<-a67
  k1[6,7]<-b67
  k0[4,5]<-a45
  k1[4,5]<-b45
  k0[5,7]<-a57
  k1[5,7]<-b57
  k012<-0.00105
  k112<-0.05596
  k12<-0.07
  intalpha<-matrix(1,nrow=8,ncol=9)
  intalpha[1,2]<-one.two
  intalpha[2,4]<-two.four
  intalpha[4,5]<-four.five
  intalpha[1,3]<-one.three
  intalpha[3,4]<-three.four
  intalpha[3,6]<-three.six
  intalpha[5,7]<-five.seven
  intalpha[6,7]<-six.seven
  intalpha[7,8]<-seven.eight
  
  trans6to7.m<-matrix(0,year2-year,(109-29))
  trans5to7.m<-matrix(0,year2-year,(109-29))
  trans6to7.f<-matrix(0,year2-year,(109-29))
  trans5to7.f<-matrix(0,year2-year,(109-29))
  fit<-sum.year<-inc.year<-malesnum<-femalesnum<-vector("list")
  for(l in 1:(109-29)){
    for(j in 1:(year2-year)){
    trans6to7.f[j, l] <- TP.f(28 + l, year + (j - 1), intalpha, int.year, mcid, f, k0, k1, k012, k112, k12)[6, 7]
    trans5to7.f[j, l] <- TP.f(28 + l, year + (j - 1), intalpha, int.year, mcid, f, k0, k1, k012, k112, k12)[5, 7]
    trans6to7.m[j, l] <- TP.m(28 + l, year + (j - 1), intalpha, int.year, mcid, f, k0, k1, k012, k112, k12)[6, 7]
    trans5to7.m[j, l] <- TP.m(28 + l, year + (j - 1), intalpha, int.year, mcid, f, k0, k1, k012, k112, k12)[5, 7]
    }
  }
  
  for(i in 1:(year2-year)){
    fit[[i]]<-Prevnum.Total(prob1[[i]],prob2[[i]],(year+i),alpha=intalpha,i=i,
                            trans6to7.f,trans5to7.f,trans6to7.m,trans5to7.m,int.year,mcid,f,female.census,male.census,
                            k0,k1,k012,k112,k12)
    prob1[[i+1]]<-fit[[i]][1:80]
    prob2[[i+1]]<-fit[[i]][81:160]
    sum.year[[i]]<-fit[[i]][161:170]
    femalesnum[[i]]<-fit[[i]][171:250]
    malesnum[[i]]<-fit[[i]][251:330]
  }
  sum<-unlist(sum.year)
  yearvar<-as.character((year+1):year2)
  table<-data.frame(cbind((yearvar),round(matrix((sum)/1000000, nrow=year2-year, ncol=10,byrow=T),2)))
  colnames(table)<-(c("Year","Normal[1]", "A[2]","N[3]","A & N[4]","MCI & A & N[5]",
                      "MCI & N[6]","Early AD[7]","Late AD[8]","Total Incidence",
                      "Incidence from State 5"))
  colnames(table) <- gsub(".", " ", colnames(table), fixed=TRUE)
  #rownames(table)[1:(year2-year)]<-paste("Year",(year+1):year2,sep="")
  tabout<-table[-(1:(int.year-2015)),]
  femalesyear=matrix(0,nrow=80, ncol=8)
  for(i in 1:80){
    femalesyear[i,]<-femalesnum[[(int.year-2017+1)]][[i]]
  }
  malesyear=matrix(0,nrow=80, ncol=8)
  for(i in 1:79){
    malesyear[i,]<-malesnum[[(int.year-2017+1)]][[i]]
  }
  totalyear=malesyear+femalesyear
  plot<-data.frame(cbind(age=c(seq(30,109,1)),totalyear/1000000))
  
  
  totalpooled=cbind(totalyear[,1:4],totalyear[,5]+totalyear[,6],totalyear[,7:8])
  plotpool<-data.frame(cbind(age=rep(c(seq(30,109,1)),6),group=c(rep(2,80), rep(3,80),rep(4,80),
                                                                 rep(5,80), rep(6,80), rep(7,80)),c(totalpooled[,2:7])/1000000))
  dpool=data.frame(lt=c(rep("solid",5),rep( "dashed", 1)))
  plotyear<-ggplot(plotpool, aes(age,V3+0.01, colour=factor(group),linetype=factor(group))) +
    geom_smooth(method="gam", fill=NA, formula=y~s(x), size=0.92)+
    #scale_y_log10(breaks =c(50,100,200,400,1000),
    #            labels = c(50,100,200,400,1000))+
    xlab("Age")+
    ylab(paste("U.S. Prevalence (millions) in",int.year))+
    scale_x_continuous(breaks = c(seq(30,110,5)), labels=c(30, "",40,"",50,"",60,"",70,"",80,"",90,"",100,"",110))+
    scale_y_continuous(breaks = c(seq(0,1.2,0.2)))+
    scale_color_manual(name="",values=c("#33CC33","#0000CC","#3399CC","#000000","#FF66CC","#CC0033"),
                       labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
    scale_linetype_manual(values=c(rep("solid",3),rep( "dotted", 1),rep("solid",2)),name="",
                          labels=c("Amyloidosis","Neurodegeneration","Amyloidosis+Neurodegen.", "MCI due to AD","Early AD", "Late AD"))+
    theme(legend.key = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key.size = unit(1.5, 'lines'),
          legend.text=element_text(size=13),text = element_text(size=12),axis.text=element_text(color = "black", size=12))
  
  return(list(tabout, plotyear))
}


original.noint <- prevnum.table(2018, 2030, 1,1,1,1,1,1,1,1,1, 1.65, 1, female.census, male.census,
                                k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                k012, k112, k12)

original.int <- prevnum.table(2018, 2030, 1,1,1,1,1,1,0.75,1,1, 1.65, 1, female.census, male.census,
                              k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                              k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                              k012, k112, k12)

fixed.noint <- prevnum.table.fixed(2018, 2030, 1,1,1,1,1,1,1,1,1, 1.65, 1, female.census, male.census,
                                k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                k012, k112, k12)

fixed.int <- prevnum.table.fixed(2018, 2030, 1,1,1,1,1,1,0.75,1,1, 1.65, 1, female.census, male.census,
                           k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                           k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                           k012, k112, k12)


original.int.57.delay <- prevnum.table(2019, 2030, 1,1,1,1,1,1,0.75,1,1, 1.65, 1, female.census, male.census, 
                                    k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                    k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                    k012, k112, k12)

fixed.int.57.delay <- prevnum.table.fixed(2019, 2030, 1,1,1,1,1,1,0.75,1,1, 1.65, 1, female.census, male.census, 
                                       k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                       k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                       k012, k112, k12)

original.int.67.delay <- prevnum.table(2019, 2030, 1,1,1,1,1,1,1,0.25,1, 1.65, 1, female.census, male.census, 
                                    k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                    k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                    k012, k112, k12)

fixed.int.67.delay <- prevnum.table.fixed(2019, 2030, 1,1,1,1,1,1,1,0.25,1, 1.65, 1, female.census, male.census, 
                                       k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
                                       k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
                                       k012, k112, k12)

# app.test1 <- prevnum.table.app(2017, 2030, 1,1,1,1,1,1,1,1,1, 1.65, 1, female.census, male.census, 
#                     k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
#                     k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
#                     k012, k112, k12)
# 
# app.int.57.delay <- prevnum.table.app(2020, 2030, 1,1,1,1,1,1,0.75,1,1, 1.65, 1, female.census, male.census, 
#                                 k0[1,2], k0[1,3], k0[3,4], k0[2,4], k0[3,6], k0[6,7], k0[4,5], k0[5,7],
#                                 k1[1,2], k1[1,3], k1[3,4], k1[2,4], k1[3,6], k1[6,7], k1[4,5], k1[5,7],
#                                 k012, k112, k12)

importants <- list(original.noint[[1]], fixed.noint[[1]], original.int[[1]], fixed.int[[1]], 
                   original.int.57.delay[[1]], fixed.int.57.delay[[1]], original.int.67.delay[[1]], fixed.int.67.delay[[1]])
names(importants) <- c("orig.noint", "fixed.noint", "orig.int", "fixed.int",
                       "orig.int.57.delay", "fixed.int.57.delay", "orig.int.67.delay", "fixed.int.67.delay")
saveRDS(importants, "fixing.forecasts.rds")


