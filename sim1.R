## Title: Simulation study 1
## 
## Description:
## This the source code to mimic the simulation study of section 4.1 in the main paper.
## 
##
## Author: Jaehyuk Jang
##
## References: 
## [1] Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." arXiv preprint arXiv:2411.10801 (2024).
## [2] Li, Fan, Laine E. Thomas, and Fan Li. "Addressing extreme propensity scores via the overlap weights." American journal of epidemiology 188.1 (2019): 250-257.
## [3] Hainmueller, Jens. "Entropy balancing for causal effects: A multivariate reweighting method to produce balanced samples in observational studies." Political analysis 20.1 (2012): 25-46.
## [4] Imai, Kosuke, and Marc Ratkovic. "Covariate balancing propensity score." Journal of the Royal Statistical Society Series B: Statistical Methodology 76.1 (2014): 243-263.
##
##
## Last Updated: 05/31/2025

############
## Libraries
############
library(rstudioapi)
library(MASS)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('utils.R')


##############
## Simulation
##############

OVL.list = list(
  strong = c(-0.5, 0.5, -0.5, 0.5, 0.5, 0.5),
  mod = c(-1, 1, -1, 0.5, -0.5, 0.5),
  weak = c(-2, 2, -2, 1, 0, 0)
)

unit.sim = function(delta, beta) {
  n = 1000
  X = mvrnorm(n, rep(0, 5), diag(1, 5))
  X. = cbind(1,X)
  ps = 1/(1+exp(-X. %*% beta))
  Z = rbinom(n, 1, prob = ps); n.t = sum(Z); n.c = n - n.t
  gamma = c(2, 2, 2, 0, 1, -1)
  tau = 1
  Y0 = (X. %*% gamma) + rnorm(n)
  Y1 = Y0 + tau
  Y = Z*Y1 + (1-Z)*Y0
  original.data = as.data.frame(cbind(Y, Z, X))
  colnames(original.data) = c("Y", "Z", "x1", "x2", "x3", "x4", "x5")
  
  # IPW
  IPW.EE.fit = ipw.sandwich_variance_estimator(X, Y, Z, weight = TRUE)
  ipw = IPW.EE.fit$est
  V.ipw = IPW.EE.fit$V %>% sqrt
  ipw_wt = IPW.EE.fit$w[Z==0]/sum(IPW.EE.fit$w[Z==0])
  
  # OW
  OW.EE.fit = ow.sandwich_variance_estimator(X, Y, Z)
  ow = OW.EE.fit$est
  V.ow = OW.EE.fit$V %>% sqrt
  
  # MIPW
  MIPW.EE.fit = mipw.sandwich_variance_estimator(X, Y, Z, delta, weight = TRUE)
  mipw = MIPW.EE.fit$est
  V.mipw = MIPW.EE.fit$V %>% sqrt
  mipw_wt = MIPW.EE.fit$w[Z==0]/sum(MIPW.EE.fit$w[Z==0])
  
  est <- c(ipw, ow, mipw, V.ipw, V.ow, V.mipw)
  
  return(est)
}

simulate = function(R, delta, beta) {
  result = matrix(NA, nrow = 6, ncol = R, dimnames = list(c("IPW", "OW", "MIPW", "V.IPW", "V.OW", "V.MIPW"), 1:R))
  for (i in 1:R) {
    sim.res <- unit.sim(delta, beta)
    result[,i] <- sim.res
    if(i%%(R/10) == 0) cat("...",i)
  }
  return(result)
}

############
## Results
############

R = 300
delta.seq = seq(0.05, 0.95, by = 0.05)
n.delta = length(delta.seq)


## Strong OVL
strong.res.array = array(NA, dim = c(6, R, n.delta))
for(i in 1:n.delta) {
  beta = OVL.list[[1]]
  d = delta.seq[i]
  strong.res.array[,,i] = simulate(R, d, beta)
  cat(paste("Completed delta:", d, "\n"))
}

## Moderate OVL
mod.res.array = array(NA, dim = c(6, R, n.delta))
for(i in 1:n.delta) {
  beta = OVL.list[[2]]
  d = delta.seq[i]
  mod.res.array[,,i] = simulate(R, d, beta)
  cat(paste("Completed delta:", d, "\n"))
}

## Weak OVL
weak.res.array = array(NA, dim = c(6, R, n.delta))
for(i in 1:n.delta) {
  beta = OVL.list[[3]]
  d = delta.seq[i]
  weak.res.array[,,i] = simulate(R, d, beta)
  cat(paste("Completed delta:", d, "\n"))
}

