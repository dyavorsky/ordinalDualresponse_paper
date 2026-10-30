# run_fits.R -- every estimation run of the empirical application, cached.
#
#   Rscript R/empirical/run_fits.R <stage> [segment]
#
# stages (each skips work whose output file already exists):
#   ingest      build the derived panel data from the SPSS export
#   aggregate   aggregate ML, clustered SEs, spec-test battery, both links
#   agglmd      aggregate marginal likelihoods by bridge sampling, both links
#   battery     the hierarchical fits: ordinal (two seeds), logit link,
#               respondent-specific cut points, choice-only, and the binary
#               fits at every cut (top-two-box and top-box in both links;
#               k = 2, 3 in logit for the report's stepwise thresholds)
#   ext         common lambda and delta in the hierarchical model
#   exthet      the same inside the respondent-specific-cut-point model
#   boot        parametric bootstrap of both lambda tests under lambda = 1,
#               simulating from the fitted common-cut model
#   boot_hetcut the same, simulating from and refitting the het-cut model
#   all         everything above, in order
#
# Fits are written to `empirical application/derived/fits/` (gitignored).
# R/empirical/analyze_cobrand.R assembles them into R/output/empirical_cobrand.rds.

source("R/lib/dgp.R"); source("R/lib/likelihood.R"); source("R/lib/compare.R")
source("R/homogeneous/spec_tests.R"); source("R/homogeneous/bayes_agg.R")
source("R/hierarchical/dgp_panel.R"); source("R/hierarchical/hb.R")
source("empirical application/ingest.R"); source("R/empirical/cluster_se.R")

SAV  <- "empirical application/Final Data 9.6.23.sav"
DER  <- "empirical application/derived"
FD   <- file.path(DER, "fits")
MCMC <- list(R = 30000, burn = 12000, thin = 10)
SEGS <- c("member", "nonmember")
dir.create(FD, showWarnings = FALSE, recursive = TRUE)

load_seg <- function(seg) readRDS(file.path(DER, sprintf("prospects_%s.rds", seg)))
fit_path <- function(name, seg) file.path(FD, sprintf("%s_%s.rds", name, seg))
cached <- function(name, seg, expr) {
  f <- fit_path(name, seg)
  if (file.exists(f)) { cat(sprintf("%-10s %-20s cached\n", seg, name)); return(invisible(readRDS(f))) }
  x <- expr
  saveRDS(x, f)
  invisible(x)
}
hb <- function(dat, seed, ...) fit_dual_hb(dat, mcmc = MCMC, seed = seed, verbose = FALSE, ...)
report_hb <- function(seg, name, fit) cat(sprintf(
  "%-10s %-20s %.1f min  accept phi %.3f zeta %s theta %s\n", seg, name, fit$runtime_min, fit$accept$phi,
  ifelse(is.na(fit$accept$zeta), "NA", sprintf("%.3f", fit$accept$zeta)),
  ifelse(is.na(fit$accept$theta), "NA", sprintf("%.3f", fit$accept$theta))))

stage_ingest <- function() {
  if (all(file.exists(file.path(DER, sprintf("prospects_%s.rds", SEGS))))) { cat("ingest cached\n"); return() }
  out <- ingest_all(SAV, DER)
  cat("label validation:", unlist(attr(out, "validation")), "\n")
}

stage_aggregate <- function() {
  f <- file.path(FD, "aggregate.rds")
  if (file.exists(f)) { cat("aggregate cached\n"); return() }
  res <- list()
  for (seg in SEGS) {
    o <- load_seg(seg); e <- o$est; el <- e; el$link <- "logit"
    mle <- fit_dual_mle(e); cl <- cluster_se_mle(mle, e); mle_l <- fit_dual_mle(el)
    st <- spec_tests_clustered(e, mle)
    res[[seg]] <- list(mle = mle, cl = cl, mle_logit = mle_l, spec = st, colnames = o$meta$colnames)
    cat(sprintf("%-10s aggregate: nll cloglog %.1f logit %.1f\n", seg, mle$nll, mle_l$nll))
  }
  saveRDS(res, f)
}

stage_agglmd <- function() {
  f <- file.path(FD, "aggregate_lmd.rds")
  if (file.exists(f)) { cat("agglmd cached\n"); return() }
  res <- list()
  for (seg in SEGS) for (lk in c("cloglog", "logit")) {
    d <- load_seg(seg)$est; d$link <- lk
    fa <- fit_dual_bayes_agg(d, iter = 20000, burn = 5000, thin = 5, seed = 7)
    res[[paste(seg, lk)]] <- c(bridge = bridge_lmd_agg(fa), laplace = as.numeric(lmd_laplace_agg(fa)),
                               accept = fa$accept, nll_mle = fa$mle$nll)
    cat(sprintf("%-10s %-8s LMD bridge %.1f\n", seg, lk, res[[paste(seg, lk)]]["bridge"]))
  }
  saveRDS(res, f)
}

stage_battery <- function() {
  for (seg in SEGS) {
    e <- load_seg(seg)$est; el <- e; el$link <- "logit"
    runs <- list(
      ord_cloglog    = list(e,  201), ord_cloglog2 = list(e, 202), ord_logit = list(el, 203),
      hetcut_cloglog = list(e,  204, het_cut = TRUE),
      mnl            = list(e,  205, choice_only = TRUE),
      t2b_logit      = list(dichotomize(el, 4L), 101), t2b_cloglog = list(dichotomize(e, 4L), 206),
      top_logit      = list(dichotomize(el, 5L), 207), top_cloglog = list(dichotomize(e, 5L), 208),
      k2_logit       = list(dichotomize(el, 2L), 209), k3_logit    = list(dichotomize(el, 3L), 210))
    for (nm in names(runs)) {
      r <- runs[[nm]]
      fit <- cached(nm, seg, do.call(hb, c(list(r[[1]], r[[2]]), r[-(1:2)])))
      if (!is.null(fit$runtime_min)) report_hb(seg, nm, fit)
    }
  }
}

stage_ext <- function(het = FALSE) {
  for (seg in SEGS) {
    e <- load_seg(seg)$est
    for (ex in c("lambda", "delta")) {
      nm <- if (het) sprintf("ext_%s_hetcut", ex) else sprintf("ext_%s", ex)
      seed <- (if (het) 401 else 301) + (ex == "delta")
      fit <- cached(nm, seg, hb(e, seed, ext = ex, het_cut = het))
      th <- fit$draws$theta
      cat(sprintf("%-10s %-20s theta %.3f sd %.3f  95%% CI [%.3f, %.3f]  null %g\n", seg, nm,
                  mean(th), sd(th), quantile(th, .025), quantile(th, .975), fit$settings$theta_null))
    }
  }
}

# Simulate from a fitted hierarchical model with lambda = 1 exactly, with
# individual parameters drawn from the estimated population distribution,
# then re-estimate lambda at the aggregate and the hierarchical level.
stage_boot <- function(het = FALSE, reps = 3L) {
  for (seg in SEGS) {
    nm <- if (het) "boot_lambda_hetcut" else "boot_lambda"
    f <- fit_path(nm, seg)
    if (file.exists(f)) { cat(sprintf("%-10s %-20s cached\n", seg, nm)); next }
    est <- load_seg(seg)$est
    base <- readRDS(fit_path(if (het) "hetcut_cloglog" else "ord_cloglog", seg))
    P <- base$settings$P; W <- base$settings$W
    phibar <- colMeans(base$draws$phibar); Sig <- apply(base$draws$Sigma, c(2, 3), mean)
    out <- data.frame(rep = seq_len(reps), lambda_agg = NA, lambda_hb = NA, lambda_hb_sd = NA)
    for (s in seq_len(reps)) {
      Phi_s <- draw_betas(est$N, phibar, Sig, seed = 1000 + s)
      B <- Phi_s[, 1:P, drop = FALSE]
      cutv <- if (het) phi_to_cutmat(Phi_s, P, W) else colMeans(base$draws$cut)
      sim <- simulate_dual_hb(est$design, B, cutv, link = "cloglog", seed = 2000 + s)
      m <- fit_dual_mle(sim)
      out$lambda_agg[s] <- fit_dual_mle_ext(sim, "lambda", m$par, 1)$theta
      fh <- hb(sim, 3000 + s, ext = "lambda", het_cut = het)
      out$lambda_hb[s] <- mean(fh$draws$theta); out$lambda_hb_sd[s] <- sd(fh$draws$theta)
      cat(sprintf("%-10s %-20s rep %d: aggregate %.3f | hierarchical %.3f (sd %.3f)\n",
                  seg, nm, s, out$lambda_agg[s], out$lambda_hb[s], out$lambda_hb_sd[s]))
    }
    saveRDS(out, f)
  }
}

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "all"
switch(stage,
  ingest      = stage_ingest(),
  aggregate   = stage_aggregate(),
  agglmd      = stage_agglmd(),
  battery     = stage_battery(),
  ext         = stage_ext(FALSE),
  exthet      = stage_ext(TRUE),
  boot        = stage_boot(FALSE),
  boot_hetcut = stage_boot(TRUE),
  all = { stage_ingest(); stage_aggregate(); stage_agglmd(); stage_battery()
          stage_ext(FALSE); stage_ext(TRUE); stage_boot(FALSE); stage_boot(TRUE) },
  stop("unknown stage: ", stage))
cat("done:", stage, "\n")
