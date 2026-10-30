# hb.R -- hierarchical Bayesian estimation of the ordinal dual-response link
# by random-walk Metropolis within Gibbs (Rossi, Allenby & McCulloch style).
#
# Model (common cut points, het_cut = FALSE):
#   beta_i ~ MVN(beta_bar, Sigma);  common cuts c (transformed to
#   zeta = (c_1, log diffs));  per-task likelihood = MNL x cumulative link
#   (paper eq-joint).
# Heterogeneous cut points (het_cut = TRUE):
#   phi_i = (beta_i, c_{i1}, log diffs_i) ~ MVN(phi_bar, Sigma) jointly.
#
# Steps per iteration:
#   1. all respondents' phi_i updated simultaneously by independent RW-MH
#      proposals with covariance s_i^2 * Sigma (vectorized likelihood).
#   2. common-cut block zeta updated by RW-MH (skipped when het_cut).
#   3. (phi_bar, Sigma) by conjugate Normal / inverse-Wishart Gibbs draws.
# Step sizes adapt toward target acceptance during burn-in only.
#
# Requires: dgp.R (panel structures), likelihood.R (ord_prob_z, par_to_cut,
# cut_to_par, row_max).

# --- vectorized respondent log-likelihoods -----------------------------------

# Individual cut matrices from the phi parameterization (het_cut case).
phi_to_cutmat <- function(Phi, P, W) {
  c1 <- Phi[, P + 1]
  if (W == 2) return(matrix(c1, ncol = 1))
  D <- exp(Phi[, (P + 2):(P + W - 1), drop = FALSE])
  cutmat <- matrix(0, nrow(Phi), W - 1)
  cutmat[, 1] <- c1
  for (w in 2:(W - 1)) cutmat[, w] <- cutmat[, w - 1] + D[, w - 1]
  cutmat
}

# N-vector of per-respondent joint log-likelihoods.
#   Bmat: N x P coefficients; cutmat: N x (W-1) utility-scale cut points.
#   ext / theta: the paper's one-parameter specification enlargements, with a
#   COMMON theta (see spec_tests.R for the aggregate versions):
#     "lambda": second-stage predictor theta * mubar      (theory: theta = 1)
#     "delta":  predictor mubar + theta * (V_{j*} - mubar) (theory: theta = 0)
resp_loglik_all <- function(dat, Bmat, cutmat, choice_only = FALSE,
                            ext = "none", theta = NULL) {
  des <- dat$design
  n <- des$n_tasks
  J <- des$J
  Bx <- Bmat[dat$respondent_id, , drop = FALSE]
  V <- matrix(rowSums(des$X * Bx), nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  mubar <- m + log(rowSums(exp(V - m)))
  Vstar <- V[cbind(seq_len(n), dat$jstar)]
  ll <- Vstar - mubar
  if (!choice_only) {
    pred <- switch(ext,
                   lambda = theta * mubar,
                   delta  = mubar + theta * (Vstar - mubar),
                   mubar)
    cut_aug <- cbind(-Inf, cutmat, Inf)
    lo <- cut_aug[cbind(dat$resp, dat$y)] - pred
    hi <- cut_aug[cbind(dat$resp, dat$y + 1L)] - pred
    ll <- ll + log(pmax(ord_prob_z(lo, hi, dat$link), 1e-312))
  }
  as.numeric(rowsum(ll, dat$resp, reorder = TRUE))
}

# --- prior and hyperparameter draws ------------------------------------------

default_hb_priors <- function(dim_phi) {
  list(
    phibar0 = rep(0, dim_phi),   # prior mean of phi_bar
    A = 1 / 100,                 # prior precision (scalar) of phi_bar
    nu = dim_phi + 3,            # IW df
    V0 = (dim_phi + 3) * diag(dim_phi),  # IW scale
    zeta_prec = 1 / 100          # prior precision of common-cut block
  )
}

# phi_bar | Sigma, Phi  ~  N(m, Vb)
draw_phibar <- function(Phi, Sigma, priors) {
  N <- nrow(Phi)
  Sinv <- chol2inv(chol(Sigma))
  prec <- N * Sinv + priors$A * diag(ncol(Phi))
  Vb <- chol2inv(chol(prec))
  m <- Vb %*% (Sinv %*% (N * colMeans(Phi)) + priors$A * priors$phibar0)
  as.numeric(m + t(chol(Vb)) %*% rnorm(ncol(Phi)))
}

# Sigma | phi_bar, Phi  ~  IW(nu + N, V0 + S)
draw_sigma <- function(Phi, phibar, priors) {
  R <- sweep(Phi, 2, phibar)
  Vn <- priors$V0 + crossprod(R)
  nun <- priors$nu + nrow(Phi)
  # inverse-Wishart draw via Wishart on the inverse scale
  Winv <- chol2inv(chol(Vn))
  Wdraw <- rWishart(1, df = nun, Sigma = Winv)[, , 1]
  chol2inv(chol(Wdraw))
}

# --- main sampler ------------------------------------------------------------

fit_dual_hb <- function(dat,
                        mcmc = list(R = 30000, burn = 10000, thin = 10),
                        het_cut = FALSE, choice_only = FALSE,
                        priors = NULL, keep_phi = TRUE,
                        ext = c("none", "lambda", "delta"),
                        seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  # Optional common extension parameter (the specification tests at the
  # hierarchical level). Starts at, and is shrunk toward, the theory's value.
  ext <- match.arg(ext)
  if (choice_only) ext <- "none"
  theta_null <- if (ext == "lambda") 1 else 0
  theta <- theta_null
  des <- dat$design
  N <- dat$N
  P <- des$P
  W <- dat$W
  dim_phi <- if (het_cut) P + W - 1 else P
  if (is.null(priors)) priors <- default_hb_priors(dim_phi)

  R_iter <- mcmc$R
  burn <- mcmc$burn
  thin <- mcmc$thin
  nkeep <- floor((R_iter - burn) / thin)

  # --- initial values
  Phi <- matrix(0, N, dim_phi)
  freq <- tabulate(dat$y, nbins = W)
  cumq <- pmin(pmax(cumsum(freq)[1:(W - 1)] / sum(freq), 1e-4), 1 - 1e-4)
  ginv <- if (dat$link == "logit") qlogis else function(q) -log(-log(q))
  cut0 <- log(des$J) + ginv(cumq)
  if (W > 2) for (k in 2:(W - 1)) cut0[k] <- max(cut0[k], cut0[k - 1] + 1e-3)
  zeta <- cut_to_par(cut0)                     # common-cut block (transformed)
  if (het_cut) Phi[, (P + 1):(P + W - 1)] <- matrix(zeta, N, W - 1, byrow = TRUE)

  phibar <- colMeans(Phi)
  Sigma <- diag(dim_phi)

  # par_to_cut(zeta, P = 0, W) reads the common-cut block directly
  cur_cutmat <- if (het_cut) phi_to_cutmat(Phi, P, W)
                else matrix(par_to_cut(zeta, 0, W), N, W - 1, byrow = TRUE)
  cur_ll <- resp_loglik_all(dat, Phi[, 1:P, drop = FALSE], cur_cutmat, choice_only,
                            ext, theta)

  # --- MH bookkeeping
  s_phi <- rep(2.93 / sqrt(dim_phi), N)
  acc_phi <- rep(0, N)
  s_zeta <- 0.05
  acc_zeta <- 0
  s_theta <- 0.05
  acc_theta <- 0
  window <- 100

  # --- storage
  keep_betabar <- matrix(NA_real_, nkeep, dim_phi)
  keep_Sigma <- array(NA_real_, c(nkeep, dim_phi, dim_phi))
  keep_cut <- if (!het_cut && !choice_only) matrix(NA_real_, nkeep, W - 1) else NULL
  keep_theta <- if (ext != "none") numeric(nkeep) else NULL
  keep_ll <- numeric(nkeep)
  keep_ll_resp <- matrix(NA_real_, nkeep, N)
  keep_Phi <- if (keep_phi) array(NA_real_, c(nkeep, N, dim_phi)) else NULL

  t0 <- Sys.time()
  ki <- 0
  for (it in seq_len(R_iter)) {

    # (1) respondent-level RW-MH, all respondents at once
    U <- chol(Sigma)
    Z <- matrix(rnorm(N * dim_phi), N, dim_phi) %*% U
    Prop <- Phi + Z * s_phi
    prop_cutmat <- if (het_cut) phi_to_cutmat(Prop, P, W) else cur_cutmat
    prop_ll <- resp_loglik_all(dat, Prop[, 1:P, drop = FALSE], prop_cutmat, choice_only,
                               ext, theta)
    Sinv <- chol2inv(chol(Sigma))
    dcur <- sweep(Phi, 2, phibar)
    dprop <- sweep(Prop, 2, phibar)
    qcur <- rowSums((dcur %*% Sinv) * dcur)
    qprop <- rowSums((dprop %*% Sinv) * dprop)
    log_alpha <- (prop_ll - 0.5 * qprop) - (cur_ll - 0.5 * qcur)
    accept <- log(runif(N)) < log_alpha
    Phi[accept, ] <- Prop[accept, ]
    cur_ll[accept] <- prop_ll[accept]
    if (het_cut && any(accept)) cur_cutmat[accept, ] <- prop_cutmat[accept, , drop = FALSE]
    acc_phi <- acc_phi + accept

    # (2) common-cut block RW-MH
    if (!het_cut && !choice_only) {
      zeta_prop <- zeta + rnorm(W - 1, sd = s_zeta)
      cut_prop <- matrix(par_to_cut(zeta_prop, 0, W), N, W - 1, byrow = TRUE)
      ll_prop <- resp_loglik_all(dat, Phi[, 1:P, drop = FALSE], cut_prop, choice_only,
                                 ext, theta)
      lp <- (sum(ll_prop) - 0.5 * priors$zeta_prec * sum(zeta_prop^2)) -
            (sum(cur_ll) - 0.5 * priors$zeta_prec * sum(zeta^2))
      if (log(runif(1)) < lp) {
        zeta <- zeta_prop
        cur_cutmat <- cut_prop
        cur_ll <- ll_prop
        acc_zeta <- acc_zeta + 1
      }
    }

    # (2b) common extension parameter by RW-MH, prior N(theta_null, 1/zeta_prec)
    if (ext != "none") {
      theta_prop <- theta + rnorm(1, sd = s_theta)
      ll_prop <- resp_loglik_all(dat, Phi[, 1:P, drop = FALSE], cur_cutmat, choice_only,
                                 ext, theta_prop)
      lp <- (sum(ll_prop) - 0.5 * priors$zeta_prec * (theta_prop - theta_null)^2) -
            (sum(cur_ll) - 0.5 * priors$zeta_prec * (theta - theta_null)^2)
      if (log(runif(1)) < lp) {
        theta <- theta_prop
        cur_ll <- ll_prop
        acc_theta <- acc_theta + 1
      }
    }

    # (3) hyperparameters
    phibar <- draw_phibar(Phi, Sigma, priors)
    Sigma <- draw_sigma(Phi, phibar, priors)

    # adapt step sizes during burn-in
    if (it <= burn && it %% window == 0) {
      rate_phi <- acc_phi / window
      s_phi <- s_phi * exp(1.5 * (rate_phi - 0.23))
      s_phi <- pmin(pmax(s_phi, 0.01), 10)
      acc_phi <- rep(0, N)
      if (!het_cut && !choice_only) {
        rate_z <- acc_zeta / window
        s_zeta <- min(max(s_zeta * exp(1.5 * (rate_z - 0.30)), 1e-3), 2)
        acc_zeta <- 0
      }
      if (ext != "none") {
        rate_t <- acc_theta / window
        s_theta <- min(max(s_theta * exp(1.5 * (rate_t - 0.30)), 1e-3), 2)
        acc_theta <- 0
      }
    }
    if (it == burn) { acc_phi <- rep(0, N); acc_zeta <- 0; acc_theta <- 0 }

    # store
    if (it > burn && (it - burn) %% thin == 0) {
      ki <- ki + 1
      keep_betabar[ki, ] <- phibar
      keep_Sigma[ki, , ] <- Sigma
      if (!is.null(keep_cut)) keep_cut[ki, ] <- par_to_cut(zeta, 0, W)
      if (!is.null(keep_theta)) keep_theta[ki] <- theta
      keep_ll[ki] <- sum(cur_ll)
      keep_ll_resp[ki, ] <- cur_ll
      if (keep_phi) keep_Phi[ki, , ] <- Phi
    }

    if (verbose && it %% 5000 == 0) {
      cat(sprintf("  iter %d/%d (%.1f min)\n", it, R_iter,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }

  post_iters <- R_iter - burn
  structure(list(
    draws = list(phibar = keep_betabar, Sigma = keep_Sigma, cut = keep_cut,
                 Phi = keep_Phi, theta = keep_theta),
    ll = keep_ll, ll_resp = keep_ll_resp,
    accept = list(phi = mean(acc_phi / post_iters),
                  zeta = if (het_cut || choice_only) NA else acc_zeta / post_iters,
                  theta = if (ext == "none") NA else acc_theta / post_iters),
    settings = list(mcmc = mcmc, het_cut = het_cut, choice_only = choice_only,
                    ext = ext, theta_null = theta_null,
                    P = P, W = W, N = N, link = dat$link, priors = priors),
    runtime_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ), class = "dr_hb_fit")
}

# Posterior-mean individual coefficients (N x P).
hb_beta_i <- function(fit) {
  apply(fit$draws$Phi[, , 1:fit$settings$P, drop = FALSE], c(2, 3), mean)
}

summary_dr_hb <- function(fit, truth = NULL) {
  pm <- colMeans(fit$draws$phibar)
  psd <- apply(fit$draws$phibar, 2, sd)
  out <- data.frame(param = paste0("phibar", seq_along(pm)), post_mean = pm, post_sd = psd)
  if (!is.null(fit$draws$cut)) {
    out <- rbind(out, data.frame(param = paste0("c", seq_len(ncol(fit$draws$cut))),
                                 post_mean = colMeans(fit$draws$cut),
                                 post_sd = apply(fit$draws$cut, 2, sd)))
  }
  if (!is.null(truth)) {
    out$truth <- truth
    out$z <- (out$post_mean - truth) / out$post_sd
  }
  out
}
