###### Generating estimates of prevalence of MCI using AN model

# setwd("/Users/Tommy/Desktop/Tommy/School/Grad School/Research/Research Brookmeyer/Code")

library(ggplot2)
library(readxl)
library(here)

rates <- read.csv(here("U.S.mortality.rates.csv"))
jack_prev_men <- as.matrix(read_excel(here("ATN_estimates_rescaled.xlsx"))[,c(1, 2, 6, 8, 9, 4, 5, 3, 7)])
men_prev <- jack_prev_men[,-1] / 100
jack_prev_fem <- as.matrix(read_excel(here("ATN_estimates_rescaled.xlsx"))[,c(1, 10, 14, 16, 17, 12, 13, 11, 15)])
fem_prev <- jack_prev_fem[,-1] / 100

source(here("BrookFuncs.R"))

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

##### For the sensitivity analysis, we also have "low" and "high" values for transitions to MCI
##### and alzheimer's

k0.low <- k0.high <- k0
k1.low <- k1.high <- k1

k0.low[4,5] <- k0[4,5] * 0.66
k0.low[5,7] <- 0.26

k0.high[4,5] <- k0[4,5] * 1.65
k0.high[5,7] <- 0.34


old.prevalences <- Prevrate.f.multi.AN(50:95, 2014, intalpha.AN, 2014, 1.65, 1, k0, k1, k012, k112, k12)

prev.preclinical <- rowSums(old.prevalences[,1:4]) / rowSums(old.prevalences[,1:8])
prev.clinical <- 1 - prev.preclinical

mci.progression <- cbind.data.frame(c(prev.clinical, prev.preclinical),
                                    rep(50:95, 2),
                                    c(rep("Clinical", length(prev.clinical)), rep("Preclinical", length(prev.clinical))))

names(mci.progression) <- c("Prevalence", "Age", "MCI Status")

ggplot(data = mci.progression, mapping = aes(x = Age, y = Prevalence)) +
  geom_line(aes(color = `MCI Status`), size = 2) +
  labs(title = "Progression of MCI Prevalence under AN Model")


# write.csv(prev.preclinical, "prev.preclinical_03.15.2021.csv")
write.csv(prev.preclinical, "prev.preclinical_06.10.2021.csv")

prevs.low <- Prevrate.f.multi.AN(50:95, 2014, intalpha.AN, 2014, 1.65, 1, k0.low, k1.low, k012, k112, k12)
prevs.high <- Prevrate.f.multi.AN(50:95, 2014, intalpha.AN, 2014, 1.65, 1, k0.high, k1.high, k012, k112, k12)

preclinical.low <- rowSums(prevs.low[,1:4]) / rowSums(prevs.low[,1:8])
preclinical.high <- rowSums(prevs.high[,1:4]) / rowSums(prevs.high[,1:8])

fem_prev_u_low <- preclinical.low * fem_prev
men_prev_u_low <- preclinical.low * men_prev
fem_prev_u_high <- preclinical.high * fem_prev
men_prev_u_high <- preclinical.high * men_prev

avg_prev_u_low <- (fem_prev_u_low + men_prev_u_low) / 2
avg_prev_u_high <- (fem_prev_u_high + men_prev_u_high) / 2


