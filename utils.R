## Title: Utility functions
## 
## Description:
## The functions here facilitate the simulation studies and real data analysis.
## Users who are only interested in the source code for the mixing strategy do not need to read this file.
##
## Author: Jaehyuk Jang
## References: 
## [1] Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." arXiv preprint arXiv:2411.10801 (2024).
## [2] Lunceford, Jared K., and Marie Davidian. "Stratification and weighting via the propensity score in estimation of causal treatment effects: a comparative study." Statistics in medicine 23.19 (2004): 2937-2960.
## [3] https://github.com/kosukeimai/CBPS
##
##
## Last Updated: 05/31/2025

source('functions.R')

## Unregister parallel clusters
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


## IPW Estimation and its sandwich variance estimator
ipw.sandwich_variance_estimator = function(X, Y, Z, weight = FALSE) {
  n = nrow(X); n.t = sum(Z)
  pi.hat = mean(Z)
  X. = cbind(1,X)
  
  ps.fit = glm(Z ~ X, family = "binomial")
  beta.hat = ps.fit$coefficients
  e = propensity_score(X., beta.hat)
  d_e = propensity_score_derivative(X., beta.hat)
  weight.hat = rep(1, n); weight.hat[Z == 0] = e[Z == 0]/(1-e[Z == 0])
  mu1.hat = sum(Z*Y)/sum(Z)
  mu0.hat = sum(weight.hat*(1-Z)*Y)/sum(weight.hat*(1-Z))
  tau.hat = mu1.hat - mu0.hat
  
  psi_1 = as.vector((Z - e) / (e * (1 - e)))*d_e
  psi_2 = Z*Y - Z*mu1.hat
  psi_3 = e / (1-e) * (1-Z) * Y - e / (1-e) * (1-Z) * mu0.hat
  
  psi = tryCatch(cbind(psi_1, psi_2, psi_3), warning = function(w) {
    return(c(psi_1, psi_2, psi_3))
  })
  meat = t(psi) %*% psi/n
  
  A.1 = as.vector((e*(1-e))^.5)*X.
  A.11 = t(A.1)%*%A.1/n
  A.31 = -t(as.vector(e/(1-e)*(1-Z)*(Y-mu0.hat)))%*%X./n
  A.22 = mean(Z)
  A.33 = mean(e/(1-e)*(1-Z))
  zero.vec = rep(0, dim(A.11)[1])
  bread = rbind(cbind(A.11, zero.vec, zero.vec),
                cbind(t(zero.vec), A.22, 0),
                cbind(A.31, 0, A.33))
  bread = solve(bread)
  Sigma = bread %*% meat %*% t(bread) / n
  p = dim(Sigma)[2]
  V = t(c(rep(0, p-2), 1, -1)) %*% Sigma %*% c(rep(0, p-2), 1, -1)
  
  out = list(V = V, est = tau.hat)
  
  if (weight == TRUE) {
    out$w = weight.hat
  }
  
  return(out)
}


## OW Estimation and its sandwich variance estimator
ow.sandwich_variance_estimator = function(X, Y, Z, weight = FALSE) {
  n = nrow(X); n.t = sum(Z)
  ps.fit = glm(Z ~ X, family = "binomial")
  X. = cbind(1, X)
  beta.hat = ps.fit$coefficients
  e = propensity_score(X., beta.hat)
  d_e = propensity_score_derivative(X., beta.hat)
  weight.hat = rep(1, n)
  weight.hat[Z == 1] = (1- e[Z == 1])/sum(1- e[Z == 1])
  weight.hat[Z == 0] = e[Z == 0]/sum(e[Z == 0])
  mu1.hat = sum(weight.hat*Z*Y)
  mu0.hat = sum(weight.hat*(1-Z)*Y)
  tau.hat = mu1.hat - mu0.hat
  
  psi_1 = as.vector((Z - e) / (e * (1 - e)))*d_e
  psi_2 = (1-e) * Z*Y - (1-e) * Z*mu1.hat
  psi_3 = e * (1-Z) * Y - e * (1-Z) * mu0.hat
  
  psi = tryCatch(cbind(psi_1, psi_2, psi_3), warning = function(w) {
    return(c(psi_1, psi_2, psi_3))
  })
  meat = t(psi) %*% psi/n
  
  
  A.1 = as.vector((e*(1-e))^.5)*X.
  A.11 = t(A.1)%*%A.1/n
  A.21 = t(as.vector(e*(1-e)*Z*(Y-mu1.hat)))%*%X./n
  A.31 = -t(as.vector(e*(1-e)*(1-Z)*(Y-mu0.hat)))%*%X./n
  A.22 = mean((1-e)*Z)
  A.33 = mean(e*(1-Z))
  zero.vec = rep(0, dim(A.11)[1])
  bread = rbind(cbind(A.11, zero.vec, zero.vec),
                cbind(A.21, A.22, 0),
                cbind(A.31, 0, A.33))
  bread = solve(bread)
  Sigma = bread %*% meat %*% t(bread) / n
  p = dim(Sigma)[2]
  V = t(c(rep(0, p-2), 1, -1)) %*% Sigma %*% c(rep(0, p-2), 1, -1)

  out = list(V = V, est = tau.hat)
  
  if (weight == TRUE) {
    out$w = weight.hat
  }
  
  return(out)
}