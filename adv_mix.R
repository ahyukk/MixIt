library(MASS)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##############
## functions
##############

## DGP
## Variation of DGP
## X1 = 0 => strong OVL
## X1 = 1 => weak OVL
dgp <- function() {
  n = 1000
  X = mvrnorm(n, rep(0, 3), diag(1, 3)) # Covariates
  X[,1] <- ifelse(X[,1] > 0, 0, 1)
  ps = 1/(1+exp(0.5 + 1.5*X[,1] - 0.5*X[,2] + 0.5*X[,3] - 1.5*X[,1]*X[,2] + 1.5*X[,1]*X[,3]))
  Z = rbinom(n, 1, prob = ps); n.t = sum(Z); n.c = n - n.t
  tau = 1                             
  Y0 = -2 + 2*(-1)^(X[,1]) + 2*X[,2] - 2*X[,3] + rnorm(n)
  Y1 = Y0 + tau
  Y = Z*Y1 + (1-Z)*Y0
  original.data = as.data.frame(cbind(Y, Z, X))
  colnames(original.data) = c("Y", "Z", "x1", "x2", "x3")
  return(original.data)
}

## e* fitting for the ATT
synth.ps.fit = function(formula, data, beta.init, delta) {
  
  mf = match.call(expand.dots = F)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf[[1L]] = as.name("model.frame")
  mf = eval(mf, parent.frame())
  mt = attr(mf, "terms")
  Z = model.response(mf, "any")
  X = if(!is.empty.model(mt)) model.matrix(mt, mf) else matrix(NA, NROW(Z), 0L)
  
  
  ## likelihood function
  ps.likelihood = function(beta.init, Z, X, delta) {
    pi.hat = sum(Z)/length(Z)
    K = (1-delta)*exp(X%*%beta.init) + delta*(pi.hat/(1-pi.hat))
    val = sum(Z*log(K) - log(1+K))
    return(val)
  }
  

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

## Simple Mixing Strategy
simp.mix <- function(data, treat, ps.formula = ps.formula, delta) {
  
  if (!is.character(treat)) {
    stop("Error: 'treat' argument must be a character.")
  }
  
  n <- dim(data)[1]
  Z <- data %>% pull(sym(treat))
  n.t <- sum(Z)
  
  mixed <- data
  
  I.t <- rbinom(n.t, 1, 1 - delta)
  
  mixed[Z==1,][I.t==0,] = data[sample(which(Z==0), sum(1 - I.t), replace = TRUE),] %>% dplyr::mutate(!!sym(treat) := 1)
  
  beta.init <- glm(ps.formula, data = data, family = "binomial")$coefficients
  synth.fit <- synth.ps.fit(ps.formula, data = mixed, beta.init = beta.init, delta = delta)
  
  mixed$ps <- synth.fit$fitted.values
  mixed$weights <- synth.fit$weights
  
  return(mixed)
  
}


## Advanced Mixing Strategy
adv.mix <- function(data, treat, var, ps.formula, d.hi, d.lo) {
  # d.hi: mix.proportion for tail group
  # d.lo: mix.proportion for body group
  # Comparison with simple mixing strategy: (d.hi + d.lo)/2 = delta
  
  if (!is.character(treat)) {
    stop("Error: 'treat' argument must be a character.")
  }
  
  n <- dim(data)[1]
  data <- separate.body_tail(data, treat, var)
  Z <- data %>% pull(sym(treat))
  tail <- data %>% pull(tail)
  tail_Z <- Z * tail
  body_Z <- Z * (1 - tail)
  
  tail_n.t <- tail_Z %>% sum
  body_n.t <- body_Z %>% sum
  
  tail_I.t <- rbinom(tail_n.t, 1, 1 - d.hi)
  body_I.t <- rbinom(body_n.t, 1, 1 - d.lo)
  
  
  mixed <- data
  
  mixed[tail_Z==1,][tail_I.t==0,] = data[sample(which(Z==0 & tail==1), sum(1 - tail_I.t), replace = TRUE),] %>% dplyr::mutate(!!sym(treat) := 1)
  mixed[body_Z==1,][body_I.t==0,] = data[sample(which(Z==0 & tail==0), sum(1 - body_I.t), replace = TRUE),] %>% dplyr::mutate(!!sym(treat) := 1)
  
  mixed.D <- mixed %>%
    dplyr::select(-c(Y, tail))
  
  beta.init <- glm(ps.formula, data = data %>% dplyr::select(-c(Y, tail)), family = "binomial")$coefficients
  
  tail_synth.fit <- synth.ps.fit(formula = ps.formula, data = mixed.D[mixed$tail == 1,], beta.init = beta.init, delta = d.hi)
  body_synth.fit <- synth.ps.fit(formula = ps.formula, data = mixed.D[mixed$tail == 0,], beta.init = beta.init, delta = d.lo)
  
  mixed$ps <- rep(1, n)
  mixed$weights <- rep(1, n)
  
  mixed$ps[mixed$tail == 1] <- tail_synth.fit$fitted.values
  mixed$ps[mixed$tail == 0] <- body_synth.fit$fitted.values
  
  mixed$weights[mixed$tail == 1] <- tail_synth.fit$weights
  mixed$weights[mixed$tail == 0] <- body_synth.fit$weights
  
  return(mixed)
}


########################
## Overlap Diagnosis
########################
set.seed(10001)

data <- dgp()

ps.hat <- glm(Z ~ x1 + x2 + x3 + x1*x2 + x1*x3, data = data, family = "binomial")$fitted.value

## Overall OVL
max.freq <- max(hist(ps.hat)$counts)
hist(ps.hat[data$Z == 1], col = rgb(0,0,1,0.2), main = "Blue: Treated, Red: Control", ylim = c(0, max.freq))
hist(ps.hat[data$Z == 0], col = rgb(1,0,0,0.2), add = T)

## X1 = 0 => strong OVL
max.freq <- max(hist(ps.hat[data$x1 == 0])$counts)
hist(ps.hat[data$Z == 1 & data$x1 == 0], col = rgb(0,0,1,0.2), main = "X1 = 0", ylim = c(0, max.freq))
hist(ps.hat[data$Z == 0 & data$x1 == 0], col = rgb(1,0,0,0.2), add = T)

## X1 = 1 => weak OVL
max.freq <- max(hist(ps.hat[data$x1 == 1])$counts)
hist(ps.hat[data$Z == 1 & data$x1 == 1], col = rgb(0,0,1,0.2), main = "X1 = 1", ylim = c(0, max.freq))
hist(ps.hat[data$Z == 0 & data$x1 == 1], col = rgb(1,0,0,0.2), add = T)
  


########################
## Estimator Comparison
########################

sim <- function(gen_data, ps.formula) {
  
  data <- gen_data()
  
  n <- dim(data)[1]
  
  ### Estimator Comparison ###
  
  ## IPW
  ps.fit <- glm(ps.formula, data = data, family = "binomial")
  ps.hat <- ps.fit$fitted.value
  w.hat <- ps.hat / (1 - ps.hat)
  EY1.1 <- mean(data$Z*data$Y)/mean(data$Z)
  ipw <- EY1.1 - mean(w.hat * (1 - data$Z) * data$Y)/mean(w.hat * (1 - data$Z))
  
  simp_w.mat <- matrix(NA, nrow = n, ncol = 200)
  adv_w.mat <- matrix(NA, nrow = n, ncol = 200)
  
  ## MIPW
  for (i in 1:200) {
    # Simple Mixing
    simp.mixed <- data %>% simp.mix(treat = "Z", ps.formula = ps.formula, delta = 0.3)
    simp_w.mat[,i] <- simp.mixed$weights
    
    # Advanced Mixing
    adv.mixed <- data %>% adv.mix(treat = "Z", var = "x1", ps.formula = ps.formula, d.hi = 0.1, d.lo = 0.5)
    adv_w.mat[,i] <- adv.mixed$weights
  }
  
  simp_w.hat <- simp_w.mat %>% apply(1, mean)
  adv_w.hat <- adv_w.mat %>% apply(1, mean)
  
  simp_mipw <- EY1.1 - mean(simp_w.hat * (1 - data$Z) * data$Y) / mean(simp_w.hat * (1 - data$Z))
  adv_mipw <- EY1.1 - mean(adv_w.hat * (1 - data$Z) * data$Y) / mean(adv_w.hat * (1 - data$Z))
  
  est.res <- c(ipw, simp_mipw, adv_mipw)
  
  ### Weight Comparison ###
  # Weight Disparity = sum of (estimated weights - estimated odds of treatment)^2
  
  ## w.pi
  w.pi <- sum(data$Z)/sum(1-data$Z)
  
  w.disparity_ipw <- sum((w.pi - w.hat[data$Z == 0])^2)                         # overall disparity of IPW
  w.disparity_ipw.0 <- sum((w.pi - w.hat[data$Z == 0 & data$x1 == 0])^2)        # x1 = 0 (strong overlap)
  w.disparity_ipw.1 <- sum((w.pi - w.hat[data$Z == 0 & data$x1 == 1])^2)        # x1 = 1 (weak overlap)
  
  w.disparity_smipw <- sum((w.pi - simp_w.hat[data$Z == 0])^2)                  # overall disparity of MIPW.simple
  w.disparity_smipw.0 <- sum((w.pi - simp_w.hat[data$Z == 0 & data$x1 == 0])^2) # x1 = 0 (strong overlap)
  w.disparity_smipw.1 <- sum((w.pi - simp_w.hat[data$Z == 0 & data$x1 == 1])^2) # x1 = 1 (weak overlap)
  
  w.disparity_amipw <- sum((w.pi - adv_w.hat[data$Z == 0])^2)                   # overall disparity of MIPW.advanced
  w.disparity_amipw.0 <- sum((w.pi - adv_w.hat[data$Z == 0 & data$x1 == 0])^2)  # x1 = 0 (strong overlap)
  w.disparity_amipw.1 <- sum((w.pi - adv_w.hat[data$Z == 0 & data$x1 == 1])^2)  # x1 = 1 (weak overlap)
  
  w.disparity = c(w.disparity_ipw, w.disparity_smipw, w.disparity_amipw)
  w.disparity.0 = c(w.disparity_ipw.0, w.disparity_smipw.0, w.disparity_amipw.0)
  w.disparity.1 = c(w.disparity_ipw.1, w.disparity_smipw.1, w.disparity_amipw.1)
  
  ### Result ###
  
  out <- c(est.res, w.disparity, w.disparity.0, w.disparity.1)
  
  return(out)
}


res <- sapply(1:1000, function(x) sim(gen_data = dgp, ps.formula = Z ~ x1 + x2 + x3 + x1*x2 + x1*x3)) %>% as.vector
sim.res <- data.frame(value = res,
                      est.type = rep(c("IPW", "MIPW.simp", "MIPW.adv"), 4),
                      data.type = c(rep("ATT", 3), rep("W.disp", 3), rep("W.disp.0", 3), rep("W.disp.1", 3)))


sapply(c("IPW", "MIPW.simp", "MIPW.adv"), function(x) {
   est <- sim.res %>% filter(data.type == "ATT" & est.type == x) %>% pull(value)
   bias <- abs(est - 1) %>% mean %>% round(3)
   se <- est %>% sd %>% round(3)
   return(c(bias, se))
 })

w.disp <- sapply(c("IPW", "MIPW.simp", "MIPW.adv"), function(x) {
  sim.res %>% filter(data.type == "W.disp" & est.type == x) %>% pull(value)
})

w.disp.0 <- sapply(c("IPW", "MIPW.simp", "MIPW.adv"), function(x) {
  sim.res %>% filter(data.type == "W.disp.0" & est.type == x) %>% pull(value) 
})

w.disp.1 <- sapply(c("IPW", "MIPW.simp", "MIPW.adv"), function(x) {
  w.disp.1 <- sim.res %>% filter(data.type == "W.disp.1" & est.type == x) %>% pull(value) 
})

boxplot(w.disp) 
boxplot(w.disp.0)
boxplot(w.disp.1)
