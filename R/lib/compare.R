# compare.R -- link-comparison primitives shared by the aggregate and
# hierarchical estimators.

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Newton-Raftery: log p(y) = -logmeanexp(-ll_draws). trim = fraction of the
# smallest-likelihood draws dropped before averaging (stabilized variant).
lmd_nr <- function(ll_draws, trim = 0) {
  if (trim > 0) {
    keep <- ll_draws >= quantile(ll_draws, trim)
    ll_draws <- ll_draws[keep]
  }
  -logmeanexp(-ll_draws)
}

# --- panel prediction helpers --------------------------------------------------
# Lifted from the notebooks (04, 08), which keep their own copies, so that the
# empirical application can source them. Requires hb.R for hb_beta_i().

# Collapse a W-point response to binary at scale point k (y >= k is "buy").
dichotomize <- function(dat, k) { dat$y <- 1L + (dat$y >= k); dat$W <- 2L; dat }

# Inclusive value of every task under respondent-specific coefficients Bmat.
task_mubar <- function(design, Bmat) {
  Bx <- Bmat[design$respondent_id, , drop = FALSE]
  V <- matrix(rowSums(design$X * Bx), ncol = design$J, byrow = TRUE)
  m <- row_max(V)
  m + log(rowSums(exp(V - m)))
}

# Share of holdout forced choices predicted by posterior-mean beta_i.
hit_rate <- function(fit, dat_holdout) {
  d <- dat_holdout$design
  Bx <- hb_beta_i(fit)[d$respondent_id, , drop = FALSE]
  V <- matrix(rowSums(d$X * Bx), ncol = d$J, byrow = TRUE)
  mean(max.col(V, ties.method = "first") == dat_holdout$jstar)
}

# Predicted Pr(y >= k) on each task of `dat`, from posterior means. For a
# binary fit (W = 2) the single cut is used whatever k is.
pred_exceed <- function(fit, dat, k) {
  mubar <- task_mubar(dat$design, hb_beta_i(fit))
  cutvec <- colMeans(fit$draws$cut)
  cc <- cutvec[min(k - 1, length(cutvec))]
  G <- if (fit$settings$link == "logit") plogis else function(z) exp(-exp(-z))
  1 - G(cc - mubar)
}
