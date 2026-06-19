library(nleqslv)
# Generate target data
dt.gen.target <- function(n_target, eta1, beta11, beta10, mu_target, p_target) {
  eta <- c(-0.1, eta1)
  beta1 <- c(0, beta11)
  beta0 <- c(0, beta10)

  intercept <- rep(1, times = n_target)
  x1 <- rsn(n_target, xi = mu_target[1], omega = 1, alpha = 0)
  x2 <- rsn(n_target, xi = mu_target[2], omega = 1, alpha = 0)
  x3 <- rbinom(n_target, size = 1, prob = p_target[1])
  x4 <- rbinom(n_target, size = 1, prob = p_target[2])
  x5 <- exp(x1 / 2)
  x6 <- exp(x2 / 2)

  dt_target <- cbind(intercept, x1, x2, x3, x4, x5, x6)
  prob_trt <- 1 / (1 + exp(-dt_target %*% eta))
  trt <- rbinom(n_target, size = 1, prob = prob_trt)

  prob_y1 <- 1 / (1 + exp(-dt_target %*% beta1))
  prob_y0 <- 1 / (1 + exp(-dt_target %*% beta0))
  y1 <- rbinom(n_target, size = 1, prob = prob_y1)
  y0 <- rbinom(n_target, size = 1, prob = prob_y0)
  y <- ifelse(trt == 1, y1, y0)

  dt_target <- cbind(y, trt, y1, y0, dt_target)
  return(dt_target)
}

# Generate source data
dt.gen.source <- function(n_source, eta1, beta11, beta10, mean_vec, p_vec, q_ind, skew_vec) {
  eta <- c(-0.1, eta1)
  beta1 <- c(0, beta11)
  beta0 <- c(0, beta10)

  intercept <- rep(1, times = n_source)
  x1 <- rsn(n_source, xi = mean_vec[q_ind, 1], omega = 1, alpha = skew_vec[q_ind, 1])
  x2 <- rsn(n_source, xi = mean_vec[q_ind, 2], omega = 1, alpha = skew_vec[q_ind, 2])
  x3 <- rbinom(n_source, size = 1, prob = p_vec[q_ind, 1])
  x4 <- rbinom(n_source, size = 1, prob = p_vec[q_ind, 2])
  x5 <- exp(x1 / 2)
  x6 <- exp(x2 / 2)

  dt_source <- cbind(intercept, x1, x2, x3, x4, x5, x6)
  prob_trt <- 1 / (1 + exp(-dt_source %*% eta))
  trt <- rbinom(n_source, size = 1, prob = prob_trt)

  prob_y1 <- 1 / (1 + exp(-dt_source %*% beta1))
  prob_y0 <- 1 / (1 + exp(-dt_source %*% beta0))
  y1 <- rbinom(n_source, size = 1, prob = prob_y1)
  y0 <- rbinom(n_source, size = 1, prob = prob_y0)
  y <- y1 * trt + y0 * (1 - trt)

  dt_source <- cbind(y, trt, y1, y0, dt_source)
  return(dt_source)
}


solve_ee <- function(ee_fn, beta_init, jac_fn = NULL) {
  # Primary: nleqslv with Newton + double-dogleg globalization + analytic Jacobian
  sol <- tryCatch(
    nleqslv(x = beta_init, fn = ee_fn, jac = jac_fn,
            method = "Newton", global = "dbldog",
            control = list(maxit = 200, ftol = 1e-8)),
    error = function(e) NULL
  )

  if (!is.null(sol) && sol$termcd %in% c(1, 2)) {
    return(sol$x)
  }

  # Fallback 1: nleqslv with Broyden (does not use jac_fn)
  sol2 <- tryCatch(
    nleqslv(x = beta_init, fn = ee_fn,
            method = "Broyden", global = "dbldog",
            control = list(maxit = 300, ftol = 1e-8)),
    error = function(e) NULL
  )

  if (!is.null(sol2) && sol2$termcd %in% c(1, 2)) {
    return(sol2$x)
  }

  # Fallback 2: try from a zero starting point (different basin)
  sol3 <- tryCatch(
    nleqslv(x = rep(0, length(beta_init)), fn = ee_fn, jac = jac_fn,
            method = "Newton", global = "dbldog",
            control = list(maxit = 200, ftol = 1e-8)),
    error = function(e) NULL
  )

  if (!is.null(sol3) && sol3$termcd %in% c(1, 2)) {
    return(sol3$x)
  }

  # Fallback 3: minimize ||ee||^2 via optim (Nelder-Mead, no gradient needed)
  obj_fn <- function(beta) sum(ee_fn(beta)^2)
  sol4 <- tryCatch(
    optim(par = beta_init, fn = obj_fn, method = "Nelder-Mead",
          control = list(maxit = 2000)),
    error = function(e) NULL
  )

  if (!is.null(sol4)) {
    return(sol4$par)
  }

  return(rep(NA_real_, length(beta_init)))
}


# For n observations with treatment vector trt and modifier vector md,
# returns a (4 x n) matrix where column i = (1, trt_i, md_i, trt_i*md_i)
build_X <- function(trt, md) {
  rbind(rep(1, length(trt)), trt, md, trt * md)
}

# --------------------------------------------------------------------------
# Target: DR estimator beta_T
# --------------------------------------------------------------------------

dr.est.target <- function(prob_fit, trt_ind, main_outcome,
                          pseudo_outcome_obs,
                          pseudo_outcome_comb,
                          pseudo_trt, modifier) {
  trt_a <- trt_ind
  yy <- main_outcome
  n <- length(yy)
  iptw_weight <- ifelse(trt_a == 1, 1 / prob_fit, 1 / (1 - prob_fit))
  iptw_weight <- pmin(iptw_weight, MAX_IPTW)     # truncate extreme weights
  md <- modifier
  pseudo_md <- rep(md, 2)

  # Compute ee1 (vectorized): P_1 augmentation term, normalized by n
  X_obs <- build_X(trt_a, md)                   # 4 x n
  resid_obs <- yy - pseudo_outcome_obs           # length n
  ee1 <- as.vector(X_obs %*% (iptw_weight * resid_obs)) / n

  # Estimating equation as function of beta (normalized by n)
  X_pseudo <- build_X(pseudo_trt, pseudo_md)     # 4 x 2n
  X_pseudo_t <- t(X_pseudo)                      # 2n x 4 (precompute)

  ee2_function <- function(beta) {
    lin_pred <- as.vector(X_pseudo_t %*% beta)
    mu_pred <- 1 / (1 + exp(-lin_pred))
    resid_pseudo <- pseudo_outcome_comb - mu_pred
    as.vector(X_pseudo %*% resid_pseudo) / n + ee1
  }

  # Analytic Jacobian: J(beta) = -(1/n) * X %*% diag(mu*(1-mu)) %*% X'
  jac_function <- function(beta) {
    lin_pred <- as.vector(X_pseudo_t %*% beta)
    mu_pred <- 1 / (1 + exp(-lin_pred))
    w <- mu_pred * (1 - mu_pred)                 # 2n x 1
    -X_pseudo %*% (w * X_pseudo_t) / n           # 4 x 4
  }

  # Initial value from GLM
  fit <- glm(yy ~ trt_a + md + md:trt_a, family = binomial(link = "logit"))
  beta_initial <- as.vector(fit$coefficients)

  if (any(is.na(beta_initial)) || !fit$converged) {
    beta_initial <- rep(0, 4)
  }

  solve_ee(ee2_function, beta_initial, jac_fn = jac_function)
}

# --------------------------------------------------------------------------
# Source: DR estimator beta_s
# --------------------------------------------------------------------------

dr.est.source <- function(beta_T, pseudo_outcome_comb, pseudo_trt, ee_source, modifier) {
  n_target <- length(pseudo_trt) / 2
  md <- modifier
  pseudo_md <- rep(md, 2)

  X_pseudo <- build_X(pseudo_trt, pseudo_md)
  ee2_function <- function(beta) {
    lin_pred <- as.vector(t(X_pseudo) %*% beta)
    mu_pred <- 1 / (1 + exp(-lin_pred))
    resid_pseudo <- pseudo_outcome_comb - mu_pred
    as.vector(X_pseudo %*% resid_pseudo) / n_target + ee_source
  }

  solve_ee(ee2_function, beta_T)
}

# --------------------------------------------------------------------------
# Target + Source: doubly robust federated estimator beta_F
# --------------------------------------------------------------------------

dr.est.comb <- function(prob_fit, trt_ind, main_outcome,
                        pseudo_outcome_obs, pseudo_outcome_comb,
                        pseudo_trt, delta, ee_source, modifier) {
  trt_a <- trt_ind
  yy <- main_outcome
  n_target <- length(yy)
  iptw_weight <- ifelse(trt_a == 1, 1 / prob_fit, 1 / (1 - prob_fit))
  iptw_weight <- pmin(iptw_weight, MAX_IPTW)     # truncate extreme weights
  md <- modifier
  pseudo_md <- rep(md, 2)

  # Compute ee1: weighted combination of target and source augmentation terms
  X_obs <- build_X(trt_a, md)
  resid_obs <- yy - pseudo_outcome_obs
  ee_target <- as.vector(X_obs %*% (iptw_weight * resid_obs)) / n_target
  ee_comb <- cbind(ee_target, ee_source)
  ee1 <- as.vector(ee_comb %*% delta)

  # Estimating equation as function of beta
  X_pseudo <- build_X(pseudo_trt, pseudo_md)
  X_pseudo_t <- t(X_pseudo)

  ee2_function <- function(beta) {
    lin_pred <- as.vector(X_pseudo_t %*% beta)
    mu_pred <- 1 / (1 + exp(-lin_pred))
    resid_pseudo <- pseudo_outcome_comb - mu_pred
    as.vector(X_pseudo %*% resid_pseudo) / n_target + ee1
  }

  jac_function <- function(beta) {
    lin_pred <- as.vector(X_pseudo_t %*% beta)
    mu_pred <- 1 / (1 + exp(-lin_pred))
    w <- mu_pred * (1 - mu_pred)
    -X_pseudo %*% (w * X_pseudo_t) / n_target
  }

  # Initial value from GLM
  fit <- glm(yy ~ trt_a + md + md:trt_a, family = binomial(link = "logit"))
  beta_initial <- as.vector(fit$coefficients)

  if (any(is.na(beta_initial)) || !fit$converged) {
    beta_initial <- rep(0, 4)
  }

  solve_ee(ee2_function, beta_initial, jac_fn = jac_function)
}

# --------------------------------------------------------------------------
# Density ratio: estimating gamma
# --------------------------------------------------------------------------

density.ratio.est <- function(x_mean_target, dt_source) {
  dt_source_covariate <- dt_source[, colnames(dt_source) %in% c("intercept", "x1", "x2", "x3", "x4")]
  n_source <- nrow(dt_source_covariate)
  X_mat <- as.matrix(dt_source_covariate)

  ee_function <- function(gam) {
    mu <- X_mat %*% gam
    w <- exp(mu) / n_source
    as.vector(t(X_mat) %*% w - x_mean_target)
  }

  solve_ee(ee_function, rep(0.1, 5))
}

# --------------------------------------------------------------------------
# Calculate residuals from the source sites
# --------------------------------------------------------------------------

SS_resi <- function(trt, y, or_fit_source, propen_fit_source, wt_source_SS, modifier) {
  n_source <- length(trt)
  iptw_weight <- ifelse(trt == 1, 1 / propen_fit_source, 1 / (1 - propen_fit_source))
  iptw_weight <- pmin(iptw_weight, MAX_IPTW)     # truncate extreme weights
  md <- modifier

  X_mat <- build_X(trt, md)                         # 4 x n_source
  combined_weight <- iptw_weight * wt_source_SS * (y - or_fit_source)
  ee_SS <- as.vector(X_mat %*% combined_weight) / n_source
  return(ee_SS)
}
