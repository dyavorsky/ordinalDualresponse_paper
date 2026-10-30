# scenarios.R -- the consultancy report's card scenarios, and sign-up
# predictions for a single offered card under each estimator.
#
# A scenario is one card configuration x shown on its own, so the inclusive
# value is S = exp(V(x)) and the probability that a respondent reports at or
# above scale point k is  Pr(y >= k | x) = 1 - G(c_{k-1} - V(x)).  "Likely to
# sign up" in the report is the top-two-box event, k = 4.
#
# Requires: R/lib/likelihood.R, R/hierarchical/hb.R, empirical application/ingest.R.

# --- the scenario table --------------------------------------------------------
# Level codes per attribute (att1..att9), in the ingest's coding. Members'
# att1 is dropped at coding time, so its value is irrelevant for them.
scenario_table <- function(segment = c("member", "nonmember")) {
  segment <- match.arg(segment)
  base <- c(att1 = 3L, att2 = 2L, att3 = 2L, att4 = 2L, att5 = 1L,
            att6 = 1L, att7 = 1L, att8 = 1L, att9 = 1L)       # the current offer
  s <- function(...) { v <- base; m <- c(...); v[names(m)] <- m; v }
  rows <- list(
    Current = s(),
    S1 = s(att4 = 1L),
    S2 = s(att4 = 1L, att6 = 2L),
    S3 = s(att4 = 1L, att6 = 2L, att8 = 2L),
    S4 = s(att4 = 1L, att6 = 2L, att8 = 2L, att5 = 2L),
    S5 = if (segment == "member") s(att5 = 2L) else s(att5 = 2L, att7 = 3L),
    S6 = if (segment == "member") s(att7 = 3L) else s(att1 = 2L),
    S7 = s(att2 = 1L, att3 = 1L, att4 = 1L),
    Worst = s(att1 = 1L, att2 = 1L, att3 = 1L, att4 = 1L, att5 = 2L,
              att6 = 2L, att7 = 3L, att8 = 2L, att9 = 2L))
  tab <- as.data.frame(do.call(rbind, rows))
  tab$scenario <- names(rows)
  tab$report_signup <- if (segment == "member")
    c(78.5, 77.4, 75.1, 72.6, 64.9, 70.8, 72.4, 64.9, 39.5) else
    c(73.7, 72.1, 66.5, 63.1, 54.1, 55.0, 65.9, 52.6, 13.9)
  rownames(tab) <- NULL
  tab[, c("scenario", paste0("att", 1:9), "report_signup")]
}

# Scenario rows in the estimation coding (same columns, same order as X).
scenario_X <- function(segment, colnames) {
  tab <- scenario_table(segment)
  X <- dummy_code(tab, drop = if (segment == "member") 1L else integer(0))
  stopifnot(identical(colnames(X), colnames))
  rownames(X) <- tab$scenario
  X
}

# --- predictors ----------------------------------------------------------------
link_cdf <- function(link) if (link == "logit") plogis else function(z) exp(-exp(-z))

# Evenly spaced subset of kept draws.
draw_index <- function(fit, ndraws) {
  nkeep <- dim(fit$draws$Phi)[1]
  unique(round(seq(1, nkeep, length.out = min(ndraws, nkeep))))
}

# Pr(y >= k | x) for each scenario row of Xs, averaged over respondents and
# posterior draws. Returns the posterior mean of the population share and its
# posterior sd (uncertainty of the aggregate, not across respondents).
# Works for common-cut and heterogeneous-cut fits.
signup_hb <- function(fit, Xs, k, ndraws = 400) {
  P <- fit$settings$P; W <- fit$settings$W
  G <- link_cdf(fit$settings$link)
  idx <- draw_index(fit, ndraws)
  out <- t(vapply(idx, function(r) {
    Phi_r <- fit$draws$Phi[r, , , drop = TRUE]
    if (is.null(dim(Phi_r))) Phi_r <- matrix(Phi_r, ncol = 1)
    V <- Phi_r[, 1:P, drop = FALSE] %*% t(Xs)                 # N x nscen
    cutk <- if (fit$settings$het_cut) phi_to_cutmat(Phi_r, P, W)[, k - 1]
            else rep(fit$draws$cut[r, k - 1], nrow(V))
    colMeans(1 - G(cutk - V))
  }, numeric(nrow(Xs))))
  list(mean = colMeans(out), sd = apply(out, 2, sd))
}

# Expected purchase probability E[p | x] = S/(S+1) for a single offered card
# under the cloglog model (the Beta(S,1) result), averaged the same way.
expected_p_hb <- function(fit, Xs, ndraws = 400) {
  stopifnot(fit$settings$link == "cloglog")
  P <- fit$settings$P
  idx <- draw_index(fit, ndraws)
  out <- t(vapply(idx, function(r) {
    Phi_r <- fit$draws$Phi[r, , , drop = TRUE]
    if (is.null(dim(Phi_r))) Phi_r <- matrix(Phi_r, ncol = 1)
    S <- exp(Phi_r[, 1:P, drop = FALSE] %*% t(Xs))
    colMeans(S / (S + 1))
  }, numeric(nrow(Xs))))
  list(mean = colMeans(out), sd = apply(out, 2, sd))
}

# The paper's eq-condexp: E[p | y = w, S] for cloglog alphas. Written with the
# ratio r = a_lo / a_hi factored out, so that an attractive task (large S)
# underflows r^S to zero instead of producing 0/0.
cond_expected_p <- function(S, alpha, w) {
  a_hi <- c(alpha, 1)[w]
  a_lo <- c(0, alpha)[w]
  lr <- if (a_lo > 0) log(a_lo / a_hi) else -Inf
  num <- -expm1((S + 1) * lr)          # 1 - r^(S+1)
  den <- -expm1(S * lr)                # 1 - r^S
  (S / (S + 1)) * a_hi * num / den
}

# Label worth across a respondent's actual tasks: for each task, the purchase
# probability the model attaches to each scale point, given the task's
# inclusive value. Uses posterior-mean beta_i and cut points (descriptive).
# Returns a tasks x W matrix of E[p | y = w, task].
label_worth_tasks <- function(fit, dat) {
  stopifnot(fit$settings$link == "cloglog", !fit$settings$het_cut)
  mubar <- task_mubar(dat$design, hb_beta_i(fit))
  alpha <- exp(-exp(-colMeans(fit$draws$cut)))
  S <- exp(mubar)
  sapply(seq_len(dat$W), function(w) cond_expected_p(S, alpha, w))
}

# --- the practitioner heuristics ---------------------------------------------
# A constant incidence from observed label shares; it cannot vary with x.
wts_juster <- c(0.05, 0.25, 0.50, 0.75, 0.95)
heur_incidence <- function(y, W = 5L, weights = wts_juster) {
  shares <- tabulate(y, W) / length(y)
  c(juster = sum(shares * weights), t2b = sum(shares[(W - 1):W]))
}
