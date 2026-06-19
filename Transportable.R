###### Effect modification:
###### Specification 1 (both PS and OR correctly specified) — BSE version, v2
###### Differences from Spec3_em_bse_v2.R:
######   - OR model: full (x1+x2+x3+x4 + interactions) — correctly specified
######   - PS model: full (x1+x2+x3+x4) everywhere — correctly specified

library(fGarch)
library(sn)
library(dplyr)
library(coxed)
library(rlist)
library(parallel)

source("fun_em.R")


### parameters setup (K, n_target, n_iter, B, batch_num, n_source_lo, n_source_hi)
K <- 10           # number of source sites
n_target <- 100
n_iter <- 20
B <- 250
batch_num <- 1
n_source_lo <- 400
n_source_hi <- 600

eta1_target <- c(0.7,0.7,-0.3,0.3,0,0)
beta11_target <- c(0.6,0.6,-0.3,0.6,0,0)   # unified with non-transportable scenario
beta10_target <- c(0.1,0.1,-0.55,0.1,0,0)

eta1_source <- c(0.7,0.7,-0.3,0.3,0,0)
beta11_source <- c(0.6,0.6,-0.3,0.6,0,0)   # unified with non-transportable scenario
beta10_source <- c(0.1,0.1,-0.55,0.1,0,0)
mu_target <- c(0.1,0.15)
p_target <- c(0.5,0.5)
mean_vec <- matrix(c(0.35,0.35,0.39,0.39,0.43,0.43), nrow = 3, byrow = TRUE)  
p_vec    <- matrix(c(0.57,0.57,0.60,0.60,0.63,0.63), nrow = 3, byrow = TRUE)
skew_vec <- matrix(c(0,0,0,0,0,0), nrow = 3, byrow = TRUE)

beta_iter <- lapply(1:n_iter, function(iter){
  set.seed(1000+iter+n_iter*(batch_num-1))
  ### generate target data
  dt_target <- dt.gen.target(n_target,eta1=eta1_target,beta11=beta11_target,beta10=beta10_target,mu_target,p_target)
  # Target: propensity score model (correct: x1+x2+x3+x4)
  propen_model_target <- glm(trt ~ x1 + x2 + x3 + x4,
                             family = binomial(link = "logit"),
                             data = as.data.frame(dt_target))
  propen_fit_target <- fitted(propen_model_target)
  if(any((propen_fit_target>0.99)|(propen_fit_target<0.01))) {
    ps_target_ind <- which((propen_fit_target>0.99)|(propen_fit_target<0.01))
    dt_target <- dt_target[-ps_target_ind,]
  }
  n_target_new <- nrow(dt_target)

  ### generate source data
  n_source_vec <- sample(n_source_lo:n_source_hi, K, replace = TRUE)
  n_source_vec <- sort(n_source_vec)
  # determine covariate parameters by tertile of source sample size
  source_quartile <- quantile(n_source_vec, probs = c(1/3,2/3))
  q_ind_quartile <- ifelse(n_source_vec<source_quartile[1], 1, ifelse(n_source_vec<source_quartile[2], 2, 3))

  dt_source_list <- lapply(1:K, function(x){
    n_source <- n_source_vec[x]
    q_ind <- q_ind_quartile[x]
    dt_source <- dt.gen.source(n_source,eta1=eta1_source,beta11=beta11_source,
                               beta10=beta10_source,mean_vec=mean_vec,p_vec=p_vec,
                               q_ind=q_ind,skew_vec=skew_vec)
    # Source: propensity score model (correct: x1+x2+x3+x4)
    propen_model_source <- glm(trt ~ x1 + x2 + x3 + x4,
                               family = binomial(link = "logit"),
                               data = as.data.frame(dt_source))
    propen_fit_source <- fitted(propen_model_source)
    if(any((propen_fit_source>0.99)|(propen_fit_source<0.01))) {
      ps_ind <- which((propen_fit_source>0.99)|(propen_fit_source<0.01))
      dt_source <- dt_source[-ps_ind,]}
    return(dt_source)
  })

  ## Target: create bootstrap samples
  dt_target_bootlist <- lapply(1:B, function(b){
    dt_target_boot <- dt_target[sort(sample(seq_len(nrow(dt_target)), nrow(dt_target), replace = TRUE)),]
    return(dt_target_boot)
  })
  dt_target_bootlist <- c(list(dt_target), dt_target_bootlist) # add the original dt_target

  ## Source: create bootstrap samples
  dt_source_bootlist <- lapply(1:K, function(k){
    dt_source <- dt_source_list[[k]]
    dt_source_boot <- lapply(1:B, function(b){
      dt_source_boots <- dt_source[sort(sample(seq_len(nrow(dt_source)), nrow(dt_source), replace = TRUE)), ]
      return(dt_source_boots)
    })
    dt_source_boot <- c(list(dt_source), dt_source_boot)
    return(dt_source_boot)
  })

  ### Target
  # Target mean vector
  dt_target_covariate_mat <- sapply(1:(B+1), function(b){
    dt_target <- dt_target_bootlist[[b]]
    dt_target_covariate <- dt_target[,colnames(dt_target)%in%c("intercept","x1","x2","x3","x4")]
    x_mean_target <- colMeans(dt_target_covariate)
    return(x_mean_target)
  })

  # Source: calculate augmentation terms
  resi_return_source <- lapply(1:K, function(k){
    dt_source_s <- dt_source_bootlist[[k]]
    ee_source_s <- lapply(1:(B+1), function(b){
      dt_source <- dt_source_s[[b]]
      # Source: propensity score model (correct: x1+x2+x3+x4)
      propen_model_source <- glm(trt ~ x1 + x2 + x3 + x4, family = binomial(link = "logit"), data = as.data.frame(dt_source))
      propen_fit_source <- fitted(propen_model_source)

      # Source: density ratio model
      x_mean_target <- dt_target_covariate_mat[,b]
      alpha_fitted <- density.ratio.est(x_mean_target,dt_source)
      dt_source_covariate <- dt_source[,colnames(dt_source)%in%c("intercept","x1","x2","x3","x4")]
      wt_source_SS <- as.vector(exp(dt_source_covariate%*%alpha_fitted))

      # Source: outcome regression model, weighted (correct: full x1+x2+x3+x4)
      outcome_model_source_SS <- glm(y ~ trt + x1 + x2 + x3 + x4 + x1:trt + x2:trt + x3:trt + x4:trt,
                                     family = quasibinomial(link = "logit"), weights = wt_source_SS,
                                     data = as.data.frame(dt_source))
      or_fit_source_SS <- fitted(outcome_model_source_SS)

      # SS estimator (weighted)
      ee_SS <- SS_resi(trt=dt_source[,"trt"],y=dt_source[,"y"],or_fit_source=or_fit_source_SS,
                       propen_fit_source=propen_fit_source,wt_source_SS=wt_source_SS,modifier=dt_source[,"x1"])

      # Source: outcome regression model (unweighted, correct: full x1+x2+x3+x4)
      outcome_model_source_SSnaive <- glm(y ~ trt + x1 + x2 + x3 + x4 + x1:trt + x2:trt + x3:trt + x4:trt,
                                          family = binomial(link = "logit"), data = as.data.frame(dt_source))
      or_fit_source_SSnaive <- fitted(outcome_model_source_SSnaive)

      # SS estimator (unweighted)
      wt_source_SSnaive <- rep(1,nrow(dt_source))
      ee_SSnaive <- SS_resi(trt=dt_source[,"trt"],y=dt_source[,"y"],or_fit_source=or_fit_source_SSnaive,
                            propen_fit_source=propen_fit_source,wt_source_SS=wt_source_SSnaive,modifier=dt_source[,"x1"])

      source_return <- cbind(ee_SS,ee_SSnaive,nrow(dt_source))
      colnames(source_return) <- c("ee_SS","ee_SSnaive","n_source")
      return(source_return)
    })
    return(ee_source_s)
  })

  # beta_T and beta_e: loop over bootstrap samples
  beta_bootlist <- lapply(1:(B+1), function(b){
    dt_target <- dt_target_bootlist[[b]]
    # Target: outcome regression model (correct: full x1+x2+x3+x4)
    outcome_model_target <- glm(y ~ trt + x1 + x2 + x3 + x4 + x1:trt + x2:trt + x3:trt + x4:trt,
                                family = binomial(link = "logit"), data = as.data.frame(dt_target))
    pseudo_y_target.obs <- fitted(outcome_model_target)
    newdata_target.trt1 <- data.frame(trt=1, dt_target[,!colnames(dt_target)%in%c("y","y0","y1","trt")])
    pseudo_y_target.trt1 <- predict(outcome_model_target, newdata = newdata_target.trt1, type = "response")
    newdata_target.trt0 <- data.frame(trt=0, dt_target[,!colnames(dt_target)%in%c("y","y0","y1","trt")])
    pseudo_y_target.trt0 <- predict(outcome_model_target, newdata = newdata_target.trt0, type = "response")

    pseudo_trt_target <- rep(c(1,0), each = nrow(dt_target))
    pseudo_outcome_target_comb <- c(pseudo_y_target.trt1, pseudo_y_target.trt0)

    # Target: propensity score model (correct: x1+x2+x3+x4)
    propen_model_target <- glm(trt ~ x1 + x2 + x3 + x4, family = binomial(link = "logit"), data = as.data.frame(dt_target))
    propen_fit_target <- fitted(propen_model_target)

    # Target-only estimator
    beta_TO <- dr.est.target(prob_fit = propen_fit_target, trt_ind = dt_target[,"trt"],
                             main_outcome = dt_target[,"y"], pseudo_outcome_obs = pseudo_y_target.obs,
                             pseudo_outcome_comb = pseudo_outcome_target_comb, pseudo_trt = pseudo_trt_target, modifier = dt_target[,"x1"])
    # Source: augmentation terms
    eq_SS_source_bootlist <- sapply(1:K, function(k){
      ee_source_s <- resi_return_source[[k]]
      ee_source_SS <- ee_source_s[[b]][,1]
      return(ee_source_SS)
    })

    eq_SSnaive_source_bootlist <- sapply(1:K, function(k){
      ee_source_s <- resi_return_source[[k]]
      ee_source_SSnaive <- ee_source_s[[b]][,2]
      return(ee_source_SSnaive)
    })

    # Source: sample size
    n_source_bootlist <- sapply(1:K, function(k){
      ee_source_s <- resi_return_source[[k]]
      n_source_boot <- ee_source_s[[b]][1,3]
      return(n_source_boot)
    })
    # calculate rho_vec
    N_total <- sum(c(nrow(dt_target),n_source_bootlist))
    N_sample <- c(nrow(dt_target),n_source_bootlist)
    rho_vec <- N_sample / N_total

    beta_SS <- dr.est.comb(prob_fit = propen_fit_target, trt_ind = dt_target[,"trt"],
                           main_outcome = dt_target[,"y"], pseudo_outcome_obs = pseudo_y_target.obs,
                           pseudo_outcome_comb = pseudo_outcome_target_comb, pseudo_trt = pseudo_trt_target,
                           delta = rho_vec, ee_source = eq_SS_source_bootlist, modifier = dt_target[,"x1"])
    beta_SSnaive <- dr.est.comb(prob_fit = propen_fit_target, trt_ind = dt_target[,"trt"],
                                main_outcome = dt_target[,"y"], pseudo_outcome_obs = pseudo_y_target.obs,
                                pseudo_outcome_comb = pseudo_outcome_target_comb, pseudo_trt = pseudo_trt_target,
                                delta = rho_vec, ee_source = eq_SSnaive_source_bootlist, modifier = dt_target[,"x1"])

    beta_Tande <- cbind(beta_TO,beta_SS,beta_SSnaive)
    colnames(beta_Tande) <- c("beta_TO", "beta_SS","beta_SSnaive")
    return(beta_Tande)
  })

  return(beta_bootlist)
})

