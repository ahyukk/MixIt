## Title: COPD Study
## 
## Description:
## This the source code to analyze the COPD dataset of Samsung Medical Center in section 5 of the main paper.
## 
##
## Author: Jaehyuk Jang
##
## References: 
## [1] Jang, Jaehyuk, Suehyun Kim, and Kwonsang Lee. "Improving Causal Estimation by Mixing Samples to Address Weak Overlap in Observational Studies." arXiv preprint arXiv:2411
## [2] Im, Yunjoo, et al. "Causal inference analysis of the radiologic progression in the chronic obstructive pulmonary disease." Scientific Reports 14.1 (2024): 17838.
##
##
## Last Updated: 06/13/2025


###########
## Libraries
###########
library(rstudioapi)
library(cobalt)
library(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('utils.R')

##############
## Utility
##############

# Covariate Balance Summary
bal.check <- function(x, bal.vars, un = FALSE) {
  bal.df <- x$Balance %>% t %>% as.data.frame %>% dplyr::select(bal.vars) %>% t
  diff.adj <- bal.df[,"Diff.Adj"] %>% as.numeric
  if (un == TRUE) {
    diff.adj <- bal.df[,"Diff.Un"] %>% as.numeric
  }
  max.asd <- diff.adj %>% abs %>% max
  max.asd.var <- rownames(bal.df)[which(abs(diff.adj) == max.asd)]
  unadj.ct_0.1 <- (abs(diff.adj) > 0.1) %>% sum
  unadj.ct_0.2 <- (abs(diff.adj) > 0.2) %>% sum
  ess <- x$Observations[2,1]
  return(c(max.asd.var, max.asd, unadj.ct_0.1, unadj.ct_0.2, ess))
}


##############
## Load Data
##############

# treated = fSAD < 0.1 (treated = 17 vs control = 61)
load(file = "data.RData")
smc <- data.smc %>% dplyr::select(age, bmi, smk, Emph, interval, Emph_outcome, treated) %>%
  mutate(Emph_change = Emph_outcome - Emph)  # Emph_change = change of emphysema progression over time

n <- dim(smc)[1]


##############
## Motivation
##############

## OVL Diagnosis: PS Modeling approach
ps.fit <- glm(treated ~ age + bmi + smk + Emph + interval, data = smc, family = "binomial")
smc$ps <- ps.fit$fitted.values
smc %>% ggplot(aes(x = ps, fill = treated)) + geom_histogram(alpha = 0.5) + theme_minimal()

## OVL Diagnosis: Weight dispersion approach
# 1. PS logistic modeling
ps.fit <- weightit(treated ~ age + bmi + smk + Emph + interval, data = smc, estimand = "ATT", method = "glm")
summary(ps.fit)
quantile(ps.fit$weights, 0.99)
sd(ps.fit$weights[smc$treated == FALSE]/sum(ps.fit$weights[smc$treated == FALSE]))
# Max Weight: 1.532 
# 99%: 1.122
# ESS: 36.2 
# SD of normalized weights: 0.0137

# 2. Entropy Balancing
ebal.fit <- weightit(treated ~ age + bmi + smk + Emph + interval, data = smc, estimand = "ATT", method = "ebal")
summary(ebal.fit)
quantile(ebal.fit$weights, 0.99)
sd(ebal.fit$weights[smc$treated == FALSE]/sum(ebal.fit$weights[smc$treated == FALSE]))
# Max Weight: 5.1213 
# 99%: 3.274 
# ESS: 37.89 
# SD of normalized weights: 0.0129 


################################
## PS Model Selection (Logistic)
################################

# Stepwise Variable Selection
null.fo <- treated ~ age + bmi + smk + Emph + interval
full.fo <- treated ~ (age + bmi + smk + Emph + interval)^2 + I(age^2) + I(bmi^2) + I(Emph^2) + I(interval^2)

null.fit <- glm(null.fo, data = smc, family = "binomial")
full.fit <- glm(full.fo, data = smc, family = "binomial")
forward.step <- step(null.fit, scope = full.fo, direction = "forward")
both.step <- step(null.fit, scope = full.fo, direction = "both")
back.step <- step(null.fit, scope = full.fo, direction = "back")

forward.step$formula
both.step$formula
back.step$formula

# Model candidates
mod0 <- null.fo # Null Model
mod1 <- both.step$formula # Stepwise variable selection (both, forward)
mod2 <- treated ~ age + bmi + smk + Emph + interval + smk*Emph # Single interaction term
mod3 <- treated ~ age + bmi + smk + Emph + interval + smk*(age + bmi + Emph + interval) # Numerous interaction terms

# Covariate Balance
bal.mod0 <- weightit(mod0, data = smc, method = "glm", estimand = "ATT") %>% bal.tab(int = TRUE, poly = 2, un = TRUE)
bal.mod1 <- weightit(mod1, data = smc, method = "glm", estimand = "ATT") %>% bal.tab(int = TRUE, poly = 2)
bal.mod2 <- weightit(mod2, data = smc, method = "glm", estimand = "ATT") %>% bal.tab(int = TRUE, poly = 2)
bal.mod3 <- weightit(mod3, data = smc, method = "glm", estimand = "ATT") %>% bal.tab(int = TRUE, poly = 2)

bal.vars <- bal.mod0$Balance %>% rownames
bal.vars <- bal.vars[-1]

mod.list <- list(glm.mod0 = bal.mod0,
                 glm.mod1 = bal.mod1,
                 glm.mod2 = bal.mod2,
                 glm.mod3 = bal.mod3)
bal.smc <- mod.list %>% sapply(function(m) bal.check(m, bal.vars)) %>% as.data.frame
rownames(bal.smc) <- c("Max.ASD.var", "Max.ASD", "#ASD>0.1", "#ASD>0.2", "ESS")
bal.unadjusted <- bal.check(bal.mod0, bal.vars = bal.vars, un = TRUE)
bal.smc <- cbind(bal.unadjusted, bal.smc)
bal.smc


#
# Summary: mod2 works the best in terms of covariate balancing
#

#######################
## Entropy Balancing
#######################

eb.fit0 <- weightit(mod0, data = smc, method = "ebal", estimand = "ATT")
summary(eb.fit0)
eb.mod0 <- eb.fit0 %>% bal.tab(int = TRUE, poly = 2) %>% bal.check(bal.vars)
bal.smc <- cbind(bal.smc, eb.mod0)

eb.fit2 <- weightit(mod2, data = smc, method = "ebal", estimand = "ATT")
summary(eb.fit2)
eb.mod2 <- eb.fit2 %>% bal.tab(int = TRUE, poly = 2) %>% bal.check(bal.vars)
bal.smc <- cbind(bal.smc, eb.mod2)
bal.smc


#######################
## Mixing Analysis
#######################

delta.seq <- seq(0.1, 0.9, by = 0.1)


# Estimation function for a single bootstrapped dataset
out.function <- function(boot.data, delta.seq) {
  
  
  print(mod2)
  print(mod0)
  
  Y <- boot.data$Emph_change
  Z <- boot.data$treated
  X <- boot.data %>% dplyr::select(-c("Emph_outcome", "Emph_change", "treated"))
  
  out <- suppressMessages(
    suppressWarnings(
      tryCatch({
        glm.wt2 <- weightit(mod2, data = boot.data, method = "glm", estimand = "ATT")$weights
        eb.wt2 <- weightit(mod2, data = boot.data, method = "ebal", estimand = "ATT")$weights
        glm.est2 <- sum(Z * glm.wt2 * Y) / sum(Z * glm.wt2) - sum((1 - Z) * glm.wt2 * Y) / sum((1 - Z) * glm.wt2)
        eb.est2 <- sum(Z * eb.wt2 * Y) / sum(Z * eb.wt2) - sum((1 - Z) * eb.wt2 * Y) / sum((1 - Z) * eb.wt2)
        
        glm.wt0 <- weightit(mod0, data = boot.data, method = "glm", estimand = "ATT")$weights
        eb.wt0 <- weightit(mod0, data = boot.data, method = "ebal", estimand = "ATT")$weights
        glm.est0 <- sum(Z * glm.wt0 * Y) / sum(Z * glm.wt0) - sum((1 - Z) * glm.wt0 * Y) / sum((1 - Z) * glm.wt0)
        eb.est0 <- sum(Z * eb.wt0 * Y) / sum(Z * eb.wt0) - sum((1 - Z) * eb.wt0 * Y) / sum((1 - Z) * eb.wt0)
        
        # model 2: treated ~ age + bmi + smk + Emph + interval + smk * Emph
        mipw.est2 <- delta.seq %>% sapply(function(d) {
          mipw.fit = mixit(mod2, outcome = "Emph_change", data = boot.data, delta = d, by = "Algorithm", weight = "glm", M = 200)
          m.est = mipw.fit$est
          return(m.est)
        })
        
        
        meb.est2 <- delta.seq %>% sapply(function(d) {
          meb.fit = mixit(mod2, outcome = "Emph_change", data = boot.data, delta = d, by = "Algorithm", weight = "ebal", M = 200)
          m.est = meb.fit$est
          return(m.est)
        })
        
        
        # model0: treated ~ age + bmi + smk + Emph + interval
        mipw.est0 <- delta.seq %>% sapply(function(d) {
          mipw.fit = mixit(mod0, outcome = "Emph_change", data = boot.data, delta = d, by = "Algorithm", weight = "glm", M = 200)
          m.est = mipw.fit$est
          return(m.est)
        })
        
        
        meb.est0 <- delta.seq %>% sapply(function(d) {
            meb.fit = mixit(mod0, outcome = "Emph_change", data = boot.data, delta = d, by = "Algorithm", weight = "ebal", M = 200)
            m.est = meb.fit$est
            return(m.est)
          })
          
        
        c(glm.est2, mipw.est2, eb.est2, meb.est2, glm.est0, mipw.est0, eb.est0, meb.est0)
      }, error = function(e) {return(rep(NA, 40))})
    )
  )
  
  return(out)
} 


# 1. Point Estimates
res <- out.function(smc, delta.seq)
print(res)

est.df <- data.frame(value = res,
                     model = rep(c("model2", "model0"), each = 20),
                     delta = rep(seq(0, 0.9, by = 0.1), times = 4),
                     Estimator = rep(c("IPW", rep("MIPW", 9), "EB", rep("MEB", 9)), times = 2))

# Model 2
est.plt2 <- est.df %>% filter(delta != 0 & model == "model2") %>% ggplot(aes(x = delta, y = value, color = Estimator)) +
  geom_point() + geom_line() +
  geom_hline(data = est.df %>% filter(delta == 0 & model == "model2"), aes(yintercept = value, color = Estimator), linetype = "dashed") +
  scale_color_manual(values = c(
    "IPW" = "#619CFF",
    "MIPW" = "#619CFF",
    "EB" = "#00BA38",
    "MEB" = "#00BA38"
  )) +
  xlab("Delta") +
  ylab("Point Estimates") +
  ggtitle("ATT Estimation") +
  coord_cartesian(ylim = c(-0.0175, 0)) +
  theme_minimal()
est.plt2

# Model 0
est.plt0 <- est.df %>% filter(delta != 0 & model == "model0") %>% ggplot(aes(x = delta, y = value, color = Estimator)) +
  geom_point() + geom_line() +
  geom_hline(data = est.df %>% filter(delta == 0 & model == "model0"), aes(yintercept = value, color = Estimator), linetype = "dashed") +
  scale_color_manual(values = c(
    "IPW" = "#619CFF",
    "MIPW" = "#619CFF",
    "EB" = "#00BA38",
    "MEB" = "#00BA38"
  )) +
  xlab("Delta") +
  ylab("Point Estimates") +
  ggtitle("ATT Estimation") +
  coord_cartesian(ylim = c(-0.0175, 0)) +
  theme_minimal()
est.plt0

ggpubr::ggarrange(est.plt2, est.plt0, ncol = 2, common.legend = TRUE)


# 2. Bootstrap Inference
B <- 2000
cl = makeCluster(4) # Parallel Computing
registerDoSNOW(cl)
pb = txtProgressBar(max = B, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

smc.boot.mat <- foreach(b = 1:B, .options.snow = opts, .combine = 'cbind') %dopar% {
  source('utils.R', local = TRUE)
  
  id.boot = sample(1:n, n, replace = TRUE)
  boot.data <- smc[id.boot,]
  
  out <- out.function(boot.data, delta.seq)
  
  return(out)
}

close(pb)
unregister_dopar()
stopCluster(cl)

smc.boot.mat %>% apply(1, function(x) sd(x, na.rm = TRUE))
smc.boot.mat %>% apply(1, function(x) sum(is.na(x)))

save(smc.boot.mat, file = "smc.boot.mat.RData")
load(file = "smc.boot.mat.RData")

se.df <- data.frame(value = smc.boot.mat %>% apply(1, function(x) sd(x, na.rm = TRUE)),
                    model = rep(c("model2", "model0"), each = 20),
                    delta = rep(seq(0, 0.9, by = 0.1), times = 4),
                    Estimator = rep(c("IPW", rep("MIPW", 9), "EB", rep("MEB", 9)), times = 2))

# Model 2
se.plt2 <- se.df %>% filter(delta != 0 & model == "model2") %>% ggplot(aes(x = delta, y = value, color = Estimator)) +
  geom_point() + geom_line() +
  geom_hline(data = se.df %>% filter(delta == 0 & model == "model2"), aes(yintercept = value, color = Estimator), linetype = "dashed") +
  scale_color_manual(values = c(
    "IPW" = "#619CFF",
    "MIPW" = "#619CFF",
    "EB" = "#00BA38",
    "MEB" = "#00BA38"
  )) +
  xlab("Delta") +
  ylab("Standard Error") +
  ggtitle("Bootstrap Inference (B = 2000)") +
  coord_cartesian(ylim = c(0.006, 0.01)) +
  theme_minimal()
se.plt2

# Model 0
se.plt0 <- se.df %>% filter(delta != 0 & model == "model0") %>% ggplot(aes(x = delta, y = value, color = Estimator)) +
  geom_point() + geom_line() +
  geom_hline(data = se.df %>% filter(delta == 0 & model == "model0"), aes(yintercept = value, color = Estimator), linetype = "dashed") +
  scale_color_manual(values = c(
    "IPW" = "#619CFF",
    "MIPW" = "#619CFF",
    "EB" = "#00BA38",
    "MEB" = "#00BA38"
  )) +
  xlab("Delta") +
  ylab("Standard Error") +
  ggtitle("Bootstrap Inference (B = 2000)") +
  coord_cartesian(ylim = c(0.006, 0.01)) +
  theme_minimal()
se.plt0

ggpubr::ggarrange(se.plt2, se.plt0, ncol = 2, common.legend = TRUE)

ggpubr::ggarrange(est.plt0, se.plt0, ncol = 2, common.legend = TRUE)
ggpubr::ggarrange(est.plt2, se.plt2, ncol = 2, common.legend = TRUE)

# 3. Covariate Balance
M <- 200
smc.delta.weight.list <- delta.seq %>% lapply(function(d) {
  Y <- smc$Emph_change
  Z <- smc$treated
  X <- smc %>% dplyr::select(-c("Emph_change", "Emph_outcome", "treated"))
  data.origin <- cbind(Z, X) %>% rename(treated = Z)
  n <- length(Z)
  n.t <- sum(Z)
  n.c <- n - n.t
  
  weights.list = lapply(seq_len(M), function(b) {
    
    m.data = mix("treated", data.origin, d) # Generated the same mixed dataset to compare the imbalance measures accurately
    
    weight.mat <- matrix(rep(1,4*n), ncol = 4)
    
    # MIPW -- model 2
    sps.fit <- synth.ps.fit(mod2, data = m.data, delta = d)
    weight.mat[Z == 0, 1] <- sps.fit$weights[Z == 0]
    
    # MEB -- model 2
    eb.fit <- suppressMessages(
      suppressWarnings(
        weightit(mod2, data = m.data, method = "ebal", estimand = "ATT")
      )
    )
    weight.mat[Z == 0, 2] = (eb.fit$weights[Z == 0] - d * (n.t/n.c)) / (1 - d)
    
    # MIPW -- model 0
    sps.fit <- synth.ps.fit(mod0, data = m.data, delta = d)
    weight.mat[Z == 0, 3] <- sps.fit$weights[Z == 0]
    
    # MEB -- model 0
    eb.fit <- suppressMessages(
      suppressWarnings(
        weightit(mod0, data = m.data, method = "ebal", estimand = "ATT")
      )
    )
    weight.mat[Z == 0, 4] = (eb.fit$weights[Z == 0] - d * (n.t/n.c)) / (1 - d)
    
    colnames(weight.mat) <- c("MIPW.2", "MEB.2", "MIPW.0", "MEB.0")
    
    return(weight.mat)
    
  })
  
  m.weights = Reduce(`+`, weights.list)/M
  
  return(m.weights)
})

delta.bal <- seq_along(delta.seq) %>% lapply(function(d) {
  temp.weight <- smc.delta.weight.list[[d]]
  bal.summary <- 1:4 %>% sapply(function(x) {
    wt <- rep(NA, n)
    wt[smc$treated == FALSE] <- temp.weight[smc$treated == FALSE,x]/sum(temp.weight[smc$treated == FALSE,x])
    wt[smc$treated == TRUE] <- temp.weight[smc$treated == TRUE,x]/sum(temp.weight[smc$treated == TRUE,x])
    smc %>% dplyr::select(-c("Emph_change", "Emph_outcome", "treated")) %>%
      bal.tab(treat = smc$treated, weights = wt, int = TRUE, poly = 2) %>%
      bal.check(bal.vars)
  })
  colnames(bal.summary) <- c("mipw.mod2", "meb.mod2", "mipw.mod0", "meb.mod0") %>% sapply(function(v) paste(v, delta.seq[d], sep = "_"))
  return(bal.summary)
})

bal.smc <- bal.smc %>% cbind(Reduce(cbind, delta.bal))
bal.smc


# 4. OVL Diagnosis: Weight Dispersion

# 4.1 Inverse Prob. Weighting -- mod2
glm.fit <- weightit(treated ~ age + bmi + smk + Emph + interval + smk*Emph, data = smc, estimand = "ATT", method = "glm")
summary(glm.fit)
quantile(glm.fit$weights, 0.99)
sd(glm.fit$weights[smc$treated == FALSE]/sum(glm.fit$weights[smc$treated == FALSE]))

# 4.2 Entropy Balancing -- mod2
ebal.fit <- weightit(treated ~ age + bmi + smk + Emph + interval + smk*Emph, data = smc, estimand = "ATT", method = "ebal")
summary(ebal.fit)
quantile(ebal.fit$weights, 0.99)
sd(ebal.fit$weights[smc$treated == FALSE]/sum(ebal.fit$weights[smc$treated == FALSE]))


# 4.3 Mixed Entropy Balancing (0.1) -- mod2
bal.tab(X, treat = Z, weights = smc.delta.weight.list[[1]][,2])
max(smc.delta.weight.list[[1]][Z == 0,2])
quantile(smc.delta.weight.list[[1]][Z == 0,2], 0.99)
sd(smc.delta.weight.list[[1]][Z == 0,2]/sum(smc.delta.weight.list[[1]][Z == 0,2]))
