## Title: MixIt
## 
## Description:
## The functions here implement the "simple mixing" strategy to estimate the ATT
## as introduced in the main paper and are planned to be further developed into an open-source R package
## by expanding the methodology into a broader framework.
## Users interested in collaboration are welcome to email jay.jaehyuk.jang@rutgers.edu.
##
## Author: Jaehyuk Jang
## References: 
## [1] Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." arXiv preprint arXiv:2411.10801 (2024).
## [2] Hainmueller, Jens. "Entropy balancing for causal effects: A multivariate reweighting method to produce balanced samples in observational studies." Political analysis 20.1 (2012): 25-46.
## [3] Imai, Kosuke, and Marc Ratkovic. "Covariate balancing propensity score." Journal of the Royal Statistical Society Series B: Statistical Methodology 76.1 (2014): 243-263.
## [4] https://github.com/kosukeimai/CBPS
##
##
## Last Updated: 06/02/2025



############
## Libraries
############
library(WeightIt)
library(CBPS)
library(numDeriv)
library(dplyr)
library(nleqslv)



######################
## Internal Functions
######################
## Function Lists
# 1. synth.ps.fit
# 2. propensity_score
# 3. propensity_score_derivative
# 4. synth_propensity_score
# 5. synth_propensity_score_derivative
# 6. find_beta
# 7. quiet


## Mixed Propensity Score Fitting (original propensity score model: logistic)
synth.ps.fit = function(formula, data, delta) {
    # formula = formula-class variable
    # data = mixed data frame
    # delta = mixing proportion ( < 1)

    ## The likelihood of mixed propensity score model
    ps.likelihood = function(beta.init, Z, X, delta) {
        pi.hat = sum(Z)/length(Z)
        K = (1-delta)*exp(X%*%beta.init) + delta*(pi.hat/(1-pi.hat))
        val = sum(Z*log(K) - log(1+K))
        return(val)
    }

    mf = match.call(expand.dots = F)
    m = match(c("formula", "data"), names(mf), 0L)
    mf = mf[c(1L, m)]
    mf[[1L]] = as.name("model.frame")
    mf = eval(mf, parent.frame())
    mt = attr(mf, "terms")
    Z = model.response(mf, "any")
    X = if(!is.empty.model(mt)) model.matrix(mt, mf) else matrix(NA, NROW(Z), 0L)
    beta.init = glm(Z ~ -1 + X, family = "binomial")$coefficients
    beta.star = optim(beta.init, ps.likelihood, Z = Z, X = X, delta = delta, control = list(fnscale = -1, maxit = 10000))$par
    names(beta.star) = colnames(X)
    w.pi = sum(Z)/(length(Z) - sum(Z))
    ps.star = 1/(1+(delta*w.pi + (1-delta)*exp(X %*% beta.star))^-1)
    w.star = rep(1, NROW(Z))
    w.star[Z == 0] = (ps.star[Z == 0]/(1-ps.star[Z == 0]) - delta*w.pi)/(1-delta)
    return(list(
        coefficients = beta.star,
        fitted.values = ps.star,
        weights = w.star
    ))
}

## Propensity score (logistic / mixed) and its derivative
propensity_score <- function(X, beta) {
  1/(1+exp(-X%*%beta))
}
propensity_score_derivative <- function(X, beta) {
  d = ifelse(is.null(ncol(X)), 1, ncol(X))
  e = propensity_score(X, beta)
  a = e * (1-e)
  as.vector(a)*X
}
synth_propensity_score <- function(X, beta, pi.hat ,delta) {
  w.pi = pi.hat/(1-pi.hat)
  1/(1+(delta*w.pi + (1-delta)*exp(X%*%beta))^(-1))
}
synth_propensity_score_derivative <- function(X, beta, pi.hat, delta) {
  d = ifelse(is.null(ncol(X)), 1, ncol(X))
  e = propensity_score(X, beta)
  d_e = propensity_score_derivative(X, beta)
  e.star = synth_propensity_score(X, beta, pi.hat, delta)
  a = (1-delta)*(1-e.star)^2/(1-e)^2
  as.vector(a)*d_e
}


## The MLE of coefficients of mixed propensity score (for M-Estimation-based method)
find_beta = function(X, Z, beta_init, delta) {
  sigma_psi = function(beta, X, Z, delta) {
    n = dim(X)[1]
    pi.hat = mean(Z)
    e.star = synth_propensity_score(X, beta, pi.hat, delta)
    d_e.star = synth_propensity_score_derivative(X, beta, pi.hat, delta)
    cons1 = delta*pi.hat+(1-pi.hat-delta)*Z - (1-pi.hat-delta*Z+pi.hat*delta)*e.star
    cons2 = (1-pi.hat)*e.star*(1-e.star)
    K.star = cons1/cons2
    psi.star = as.vector(K.star)*d_e.star
    sum_psi = colSums(psi.star)
    return(sum_psi)
  }
  root_result = tryCatch({nleqslv(beta_init, function(beta) sigma_psi(beta, X, Z, delta), control = list(maxit = 10000))},
                         error = function(e) stop("Non-finite value in estimating coefficients of mixed propensity score"))
  return(root_result$x)
}

## Ignore error, warning, any additional messages
quiet = function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}

######################
## Top-level Functions
######################
# 1. mix
# 2. mipw.sandwich_variance_estimator
# 3. mixit

## (Simple) Mixing Algorithm
mix = function(treat, data, delta) {
  # treat = binary variable
  # data = dataframe
  # delta = mixing proportion ( < 1)
  treat.vec = tryCatch({
    data %>% dplyr::pull(treat)
  }, error = function(e) {
    stop(paste("Treatment variable '", treat, "' is not found in the dataset"))
  })
  
  mix.data = data
  n.t = sum(treat.vec)
  n.c = length(treat.vec) - n.t
  I.delta = rbinom(n.t, 1, prob = delta)
  mix.data[treat.vec == 1, ][I.delta == 1, ] = data[treat.vec == 0, ][sample(1:n.c, sum(I.delta), replace = T), ]
  mix.data[[treat]] = treat.vec
  
  return(mix.data)
}



## M-Estimation-based Method to estimate MIPW
mipw.sandwich_variance_estimator = function(X, Y, Z, delta, weight = FALSE) {
  n = dim(X)[1]
  pi.hat = mean(Z)
  ps.fit = glm(Z ~ X, family = "binomial")
  beta.init = ps.fit$coefficients
  X. = cbind(1,X)
  
  beta.hat = find_beta(X., Z, beta.init, delta)
  e = propensity_score(X., beta.hat)
  d_e = propensity_score_derivative(X., beta.hat)
  e.star = synth_propensity_score(X., beta.hat, pi.hat, delta)
  d_e.star = synth_propensity_score_derivative(X., beta.hat, pi.hat, delta)
  tilde_e.star = delta*(1-e)^2/(1-pi.hat+pi.hat*delta-delta*e)^2
  tilde_d_e.star = as.vector(-2*(1-delta)*(1-e.star)/(1-e)^2*tilde_e.star)*d_e
  
  mu1.hat = sum(Z*Y)/sum(Z)
  weight.hat = rep(1, n); weight.hat[Z == 0] = e[Z == 0]/(1-e[Z == 0])
  mu0.hat = sum(weight.hat*(1-Z)*Y)/sum(weight.hat*(1-Z))
  tau.hat = mu1.hat - mu0.hat
  
  psi_0 = Z - pi.hat
  cons1 = delta*pi.hat+(1-pi.hat-delta)*Z - (1-pi.hat-delta*Z+pi.hat*delta)*e.star
  cons2 = (1-pi.hat)*e.star*(1-e.star)
  K.star = cons1/cons2
  psi_1 = as.vector(K.star)*d_e.star
  psi_2 = Z*Y - Z*mu1.hat
  psi_3 = e/(1-e)*(1-Z)*Y - e/(1-e)*(1-Z)*mu0.hat
  psi = cbind(psi_0, psi_1, psi_2, psi_3)
  meat = t(psi)%*%psi/n
  
  phi = (1-pi.hat)*(1-2*e.star)*tilde_e.star - e.star*(1-e.star)
  L.star = (cons2*(delta-Z + (1-delta)*e.star - (1-pi.hat-delta*Z+pi.hat*delta)*tilde_e.star) - phi*cons1)/cons2^2
  A.1 = as.vector(L.star)*d_e.star + as.vector(K.star)*tilde_d_e.star
  
  J.star = ((1-pi.hat-delta*Z+pi.hat*delta)*e.star^2+(delta*pi.hat+(1-pi.hat-delta)*Z)*(1-2*e.star))/((1-pi.hat)*e.star^2*(1-e.star)^2)
  A.2_1 = as.vector(J.star^.5)*d_e.star
  D.star= (1-delta)*(1-pi.hat)^2*(2*delta+(1-pi.hat+pi.hat*delta-delta*e)*e*(1-e)*(1-2*e))/(1-pi.hat+pi.hat*delta-delta*e)^3
  
  zero.vector = rep(0, dim(X.)[2])
  A.10 = -colMeans(A.1)
  
  A.11 = t(A.2_1)%*%A.2_1/n - t(as.vector(K.star)*d_e)%*%(as.vector(D.star)*d_e)/n
  A.31 = -t(as.vector((e/(1-e))*(1-Z)*(Y-mu0.hat)))%*%X./n
  A.22 = mean(Z)
  A.33 = mean((e/(1-e))*(1-Z))
  V = tryCatch(
    {
      bread1 = rbind(cbind(1, t(zero.vector), 0, 0),
                     cbind(A.10, A.11, zero.vector, zero.vector),
                     cbind(0, t(zero.vector), A.22, 0),
                     cbind(0, A.31, 0, A.33))
      bread = ginv(bread1)
      Sigma = bread %*% meat %*% t(bread) / n
      p = dim(Sigma)[1]
      V = t(c(rep(0, p-2), 1, -1)) %*% Sigma %*% c(rep(0, p-2), 1, -1)
    },
    error=function(e) {
      message('The Jacobian matrix is not invertible.')
      return(NA)
    }
    
  )
  
  out = list(V = V, est = tau.hat)
  
  if (weight == TRUE) {
    out$w = weight.hat
  }
  
  return(out)
}


## (simple) mixing strategy to estimate the ATT
mixit = function(formula, outcome, data, delta, by, weight = NULL, M = 200, ...) {
    # formula = formula-class variable (response ~ covariates)
    # treat = column name of the observed outcome variable in data
    # data = dataframe-type observed sample
    # delta = mixing proportion ( < 1)
    # by = implementation method (M-estimation or Algorithm)
    # weight = weighting method via Mixing Algorithm-based method (glm, ebal, cbps)
    # M = mixing algorithm resampling iterations
    # ... = additional arguments are passed to `ebal` and `cbps`
  
  mf = match.call(expand.dots = F)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf[[1L]] = as.name("model.frame")
  mf = eval(mf, parent.frame())
  mt = attr(mf, "terms")
  Z = model.response(mf, "any")
  X = if(!is.empty.model(mt)) model.matrix(mt, mf) else matrix(NA, NROW(Z), 0L)
  Y = tryCatch({
    data %>% dplyr::pull(outcome)
  }, error = function(e) {
    stop(paste("The observed outcome variable '", outcome, "' is not found in the dataset"))
  })
  
  if (by == "M-estimation") {
    
    out = mipw.sandwich_variance_estimator(X = X, Y = Y, Z = Z, delta = delta, weight = TRUE)
          
  } else if (by == "Algorithm") {
    
    original.data <- model.frame(formula, data = data)
    n <- length(Z)
    n.t <- sum(Z)
    n.c <- sum(1-Z)
    
    if (weight %in% c("glm", "ebal", "CBPS")) {
      
      if (weight == "glm") {
        
        mix.weights <- foreach(m = 1:M, .combine = `+`) %do% {
          mix.data <- mix(treat = deparse(formula[[2]]), data = original.data, delta = delta)
          sps.fit <- synth.ps.fit(formula, data = mix.data, delta = delta)
          return(sps.fit$weights[Z == 0])
        } / M
        
        
        
      } else if (weight == "ebal") {
        
        mix.weights <- foreach(m = 1:M, .combine = `+`) %do% {
          mix.data <- mix(treat = deparse(formula[[2]]), data = original.data, delta = delta)
          eb.fit <- weightit(formula, data = mix.data, method = "ebal", estimand = "ATT", ...) %>% suppressWarnings %>% suppressMessages
          return((eb.fit$weights[Z == 0] - delta * (n.t / n.c)) / (1 - delta))
        } / M
        
        
      } else if (weight == "CBPS") {
        
        mix.weights <- foreach(m = 1:M, .combine = `+`) %do% {
          mix.data <- mix(treat = deparse(formula[[2]]), data = original.data, delta = delta)
          cbps.fit <- CBPS(formula, data = mix.data, standardize = FALSE, ...) %>% quiet
          return(((n.t / n) * cbps.fit$weights[Z==0] - delta * (n.t / n.c) ) / (1 - delta))
        } / M
        
        
      }
      
      temp.weights = rep(1, n)
      temp.weights[Z == 0] = mix.weights
      mix.weights = temp.weights
      
      mix.est = sum(Z * mix.weights * Y) / sum(Z * mix.weights) - sum((1 - Z) * mix.weights * Y) / sum((1-Z) * mix.weights)
      
      out <- list(
        est = mix.est,
        weights = mix.weights
      )
      
    } else {
      
      stop("The only available applicable weighting methods are: 'glm', 'ebal' and 'CBPS'")
      
    }
    
  } else {
    
    stop("The only available mixing implementations are: 'M-estimation' and 'Algorithm'")
    
  }
  
  
  
  return(out)
  
}