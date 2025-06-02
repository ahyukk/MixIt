## Title: Simulation study 2
## 
## Description:
## This the source code to mimic the simulation study of section 4.2 in the main paper.
## 
##
## References: 
## [1] Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." arXiv preprint arXiv:2411.10801 (2024).
## [2] Kang, Joseph DY, and Joseph L. Schafer. "Demystifying double robustness: A comparison of alternative strategies for estimating a population mean from incomplete data." (2007): 523-539.
##
##
## Last Updated: 06/02/2025

############
## Libraries
############
library(rstudioapi)
library(MASS)
library(doSNOW)
library(ebal)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('utils.R')


##############
## Simulation
##############
unit.sim = function(delta, setting) {
  n = 1000
  X = mvrnorm(n, rep(0, 5), diag(1, 5))
  X. = cbind(1,X)
  ps = 1/(1+exp(-X. %*% c(-2, 2, -2, 1, 0, 0)))
  Z = rbinom(n, 1, prob = ps); n.t = sum(Z); n.c = n - n.t
  

  if (setting == 1) {
    # Correctly Specified
    gamma = c(2, 2, 2, 0, 1, -1)
    tau = 1
    Y0 = (X. %*% gamma) + rnorm(n)
    Y1 = Y0 + tau
    Y = Z*Y1 + (1-Z)*Y0
    original.data = as.data.frame(cbind(Y, Z, X))
    colnames(original.data) = c("Y", "Z", "x1", "x2", "x3", "x4", "x5")

  } else if (setting == 2) {
    # Misspecified (Kang & Schafer, 2007)
    gamma = c(-13.7, 27.4, 13.7, 13.7, 13.7, 13.7)
    tau = 210
    Y0 = (X. %*% gamma) + rnorm(n)
    Y1 = Y0  + tau
    Y = Z*Y1 + (1-Z)*Y0
    D = cbind(exp(X[,1]/2),
              X[,2]/(1+exp(X[,1])) + 10,
              (X[,1]*X[,3]/25+0.6)^3,
              (X[,1]+X[,4]+2)^2,
              abs(X[,3]-X[,5]+1)^.5)
    
    original.data = as.data.frame(cbind(Y, Z, D))
    colnames(original.data) = c("Y", "Z", "x1", "x2", "x3", "x4", "x5")
    
    
  }
  
  # EB
  eb.fit = weightit(Z ~ x1 + x2 + x3 + x4 + x5, data = original.data, method = "ebal", estimand = "ATT") %>% suppressWarnings
  eb.est = sum(Z * eb.fit$weights * Y) / sum(Z * eb.fit$weights) - sum((1 - Z) * eb.fit$weights * Y) / sum((1 - Z) * eb.fit$weights)
  
  # OW
  ow.fit = weightit(Z ~ x1 + x2 + x3 + x4 + x5, data = original.data, method = "glm", estimand = "ATO") %>% suppressWarnings
  ow.est = sum(Z * ow.fit$weights * Y) / sum(Z * ow.fit$weights) - sum((1 - Z) * ow.fit$weights * Y) / sum((1 - Z) * ow.fit$weights)
  
  # MIPW
  meb.fit = mixit(Z ~ x1 + x2 + x3 + x4 + x5, outcome = "Y", data = original.data, delta = delta, by = "Algorithm", weight = "ebal", M = 200)
  meb.est = meb.fit$est
  
  est <- c(eb.est, ow.est, meb.est)
  
  return(est)
}

simulate = function(R, delta, setting) {
  result = matrix(NA, nrow = 3, ncol = R, dimnames = list(c("EB", "OW", "MEB"), 1:R))
  for (i in 1:R) {
    sim.res <- unit.sim(delta, setting)
    result[,i] <- sim.res
    if(i%%(R/10) == 0) cat("...",i)
  }
  return(result)
}

############
## Results
############

R = 3000
delta.seq = seq(0.05, 0.95, by = 0.05)
n.delta = length(delta.seq)


## Correctly specified
correct.res.array = array(NA, dim = c(3, R, n.delta))
for(i in 1:n.delta) {
  d = delta.seq[i]
  correct.res.array[,,i] = simulate(R, d, setting = 1)
  cat(paste("Completed delta:", d, "\n"))
}

## Misspecified
incorrect.res.array = array(NA, dim = c(3, R, n.delta))
for(i in 1:n.delta) {
  d = delta.seq[i]
  incorrect.res.array[,,i] = simulate(R, d, setting = 2)
  cat(paste("Completed delta:", d, "\n"))
}
