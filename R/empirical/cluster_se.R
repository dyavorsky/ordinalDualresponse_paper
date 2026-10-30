# cluster_se.R -- respondent-clustered (sandwich) standard errors for the
# aggregate estimators on panel data.
#
# fit_dual_mle() and fit_dual_mle_ext() report inverse-Hessian standard errors,
# which treat every task as an independent observation. On real panel data the
# twelve tasks of one respondent share her tastes, so those SEs are too small
# and the Wald tests built on them reject too often. The sandwich
#   V = H^{-1} (sum_i s_i s_i') H^{-1}
# with s_i the score of respondent i's log-likelihood is valid under
# within-respondent dependence and under misspecification of the aggregate
# model, which is the situation the aggregate fits are in.
#
# Likelihood-ratio statistics are left as they are; they have no clustered
# analogue here, and the write-up says so.
#
# Requires: R/lib/likelihood.R, R/homogeneous/spec_tests.R.

# --- per-respondent negative log-likelihoods ---------------------------------

# Same algebra as negloglik(), summed within respondent instead of overall.
negloglik_by_resp <- function(par, dat) {
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  beta <- par[1:P]
  cut <- par_to_cut(par, P, W)
  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  mubar <- m + log(rowSums(exp(V - m)))
  ll <- V[cbind(seq_len(n), dat$jstar)] - mubar +
        log(pmax(ord_prob(mubar, cut, dat$y, dat$link), 1e-312))
  -as.numeric(rowsum(ll, dat$resp, reorder = TRUE))
}

# The lambda / delta extensions, per respondent. Layout follows negloglik_ext().
negloglik_ext_by_resp <- function(par, dat, ext = c("lambda", "delta"), fix = NULL) {
  ext <- match.arg(ext)
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  if (is.null(fix)) {
    theta <- par[P + W]
    par_core <- par[1:(P + W - 1)]
  } else {
    theta <- fix
    par_core <- par
  }
  beta <- par_core[1:P]
  cut <- par_to_cut(par_core, P, W)
  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  mubar <- m + log(rowSums(exp(V - m)))
  Vstar <- V[cbind(seq_len(n), dat$jstar)]
  pred <- if (ext == "lambda") theta * mubar else mubar + theta * (Vstar - mubar)
  ll <- Vstar - mubar + log(pmax(ord_prob(pred, cut, dat$y, dat$link), 1e-312))
  -as.numeric(rowsum(ll, dat$resp, reorder = TRUE))
}

# --- the sandwich -------------------------------------------------------------

# Jacobian of a vector-valued function by central differences: N x k.
num_jac <- function(fvec, x, eps = 1e-6) {
  f0 <- fvec(x)
  J <- matrix(NA_real_, length(f0), length(x))
  for (k in seq_along(x)) {
    h <- eps * max(1, abs(x[k]))
    xp <- x; xp[k] <- xp[k] + h
    xm <- x; xm[k] <- xm[k] - h
    J[, k] <- (fvec(xp) - fvec(xm)) / (2 * h)
  }
  J
}

# Clustered covariance of the parameter vector at `par`.
#   fvec: function returning per-respondent negative log-likelihoods
#   H:    observed information (Hessian of the total negative log-likelihood)
cluster_vcov <- function(par, fvec, H) {
  S <- num_jac(fvec, par)                 # respondent scores (of the nll)
  N <- nrow(S)
  meat <- crossprod(S) * N / (N - 1)      # small-sample factor
  Hinv <- solve(H)
  Hinv %*% meat %*% Hinv
}

# Clustered SEs for a fit_dual_mle() result: beta, cut points (delta method),
# and alpha = F(cut) under the cloglog link.
cluster_se_mle <- function(fit, dat) {
  P <- dat$design$P
  W <- dat$W
  fn <- function(p) negloglik(p, dat)
  H <- optimHess(fit$par, fn)
  V <- cluster_vcov(fit$par, function(p) negloglik_by_resp(p, dat), H)
  Jc <- matrix(0, W - 1, length(fit$par))
  Jc[, P + 1] <- 1
  if (W > 2) {
    d <- exp(fit$par[(P + 2):(P + W - 1)])
    for (w in 2:(W - 1)) Jc[w, (P + 2):(P + w)] <- d[1:(w - 1)]
  }
  out <- list(vcov = V,
              se_beta = sqrt(diag(V)[1:P]),
              se_cut = sqrt(diag(Jc %*% V %*% t(Jc))))
  if (dat$link == "cloglog") out$se_alpha <- out$se_cut * fit$alpha * (-log(fit$alpha))
  out
}

# The specification-test battery with clustered Wald statistics added.
# Columns from spec_test_battery() are kept; `se_cl`, `wald_z_cl`, `p_wald_cl`
# are the clustered versions. LR columns are unchanged and unclustered.
spec_tests_clustered <- function(dat, mle = NULL) {
  if (is.null(mle)) mle <- fit_dual_mle(dat)
  batt <- spec_test_battery(dat, mle)
  start_core <- mle$par
  cl_theta_se <- function(ext, theta0) {
    f <- fit_dual_mle_ext(dat, ext, start_core, theta0)
    fn <- function(p) negloglik_ext(p, dat, ext = ext)
    H <- optimHess(f$par, fn)
    V <- cluster_vcov(f$par, function(p) negloglik_ext_by_resp(p, dat, ext = ext), H)
    k <- length(f$par)
    c(theta = f$theta, se_cl = sqrt(V[k, k]))
  }
  l <- cl_theta_se("lambda", 1)
  d <- cl_theta_se("delta", 0)
  batt$se_cl <- c(l["se_cl"], d["se_cl"], d["se_cl"])
  batt$wald_z_cl <- c((l["theta"] - 1) / l["se_cl"],
                      d["theta"] / d["se_cl"],
                      (d["theta"] - 1) / d["se_cl"])
  batt$p_wald_cl <- 2 * pnorm(-abs(batt$wald_z_cl))
  rownames(batt) <- NULL
  batt
}
