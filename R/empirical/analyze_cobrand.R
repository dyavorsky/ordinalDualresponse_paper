# analyze_cobrand.R -- assemble the empirical application's results from the
# cached fits into R/output/empirical_cobrand.rds and R/output/figs/.
#
# Everything written here is blind and respondent-free: population-level
# estimates, posterior summaries, scenario shares, test statistics. The
# respondent-level data and the fit objects stay under the gitignored
# `empirical application/derived/`.
#
# Run from the project root:  Rscript R/empirical/analyze_cobrand.R
# Requires the fits produced by the battery, ext and aggregate scripts (the
# notebook documents them); a missing fit is reported and skipped.

source("R/lib/dgp.R"); source("R/lib/likelihood.R"); source("R/lib/compare.R")
source("R/homogeneous/spec_tests.R"); source("R/hierarchical/dgp_panel.R")
source("R/hierarchical/hb.R"); source("R/hierarchical/compare_hb.R")
source("empirical application/ingest.R"); source("R/empirical/cluster_se.R")
source("R/empirical/scenarios.R")

DER <- "empirical application/derived"
FD  <- file.path(DER, "fits")
OUT <- file.path("R", "output")
FIG <- file.path(OUT, "figs")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

load_fit <- function(name, seg) {
  f <- file.path(FD, sprintf("%s_%s.rds", name, seg))
  if (!file.exists(f)) { message("missing fit: ", f); return(NULL) }
  readRDS(f)
}
pm <- function(x) colMeans(x)
psd <- function(x) apply(x, 2, sd)
Fclog <- function(c) exp(-exp(-c))

agg_all <- readRDS(file.path(FD, "aggregate.rds"))
agg_lmd <- if (file.exists(file.path(FD, "aggregate_lmd.rds")))
  readRDS(file.path(FD, "aggregate_lmd.rds")) else NULL

results <- list()
for (seg in c("member", "nonmember")) {
  message("== ", seg)
  o <- readRDS(file.path(DER, sprintf("prospects_%s.rds", seg)))
  est <- o$est; hold <- o$holdout; full <- o$full; meta <- o$meta
  cn <- meta$colnames
  R <- list(segment = seg, N = meta$N, P = meta$P, colnames = cn,
            attributes = lapply(meta$attributes, function(a) a[c("var", "label", "levels")]),
            T_est = est$design$T_tasks, T_hold = hold$design$T_tasks,
            duplicates = meta$duplicates[c("n_tasks", "n_total")])

  # ---- A. descriptives --------------------------------------------------------
  tt <- rep(seq_len(full$design$T_tasks), times = full$N)
  R$descr <- list(
    y_freq = tabulate(full$y, 5L), y_share = tabulate(full$y, 5L) / length(full$y),
    mean_y_by_task = tapply(full$y, tt, mean),
    t2b_by_task = tapply(full$y >= 4L, tt, mean),
    jstar_share = tabulate(full$jstar, 3L) / length(full$jstar))

  # ---- B. aggregate ML ----------------------------------------------------------
  a <- agg_all[[seg]]
  R$agg <- list(
    table = data.frame(coef = cn, beta = a$mle$beta, se = a$mle$se_beta,
                       se_cl = a$cl$se_beta, row.names = NULL),
    cut = a$mle$cut, se_cut = a$mle$se_cut, se_cut_cl = a$cl$se_cut,
    alpha = a$mle$alpha, se_alpha = a$mle$se_alpha, se_alpha_cl = a$cl$se_alpha,
    nll_cloglog = a$mle$nll, nll_logit = a$mle_logit$nll,
    beta_logit = a$mle_logit$beta, se_beta_logit = a$mle_logit$se_beta,
    spec = a$spec,
    lmd = if (!is.null(agg_lmd)) list(cloglog = agg_lmd[[paste(seg, "cloglog")]],
                                      logit = agg_lmd[[paste(seg, "logit")]]) else NULL)

  # ---- C. hierarchical, common cut points ----------------------------------
  ord  <- load_fit("ord_cloglog", seg); ord2 <- load_fit("ord_cloglog2", seg)
  ordl <- load_fit("ord_logit", seg);   mnl  <- load_fit("mnl", seg)
  if (!is.null(ord)) {
    R$hb <- list(
      phibar = data.frame(coef = cn, mean = pm(ord$draws$phibar), sd = psd(ord$draws$phibar), row.names = NULL),
      cut = pm(ord$draws$cut), cut_sd = psd(ord$draws$cut),
      alpha = Fclog(pm(ord$draws$cut)),
      alpha_ci = apply(Fclog(ord$draws$cut), 2, quantile, c(.025, .975)),
      Sigma_diag = pm(t(apply(ord$draws$Sigma, 1, diag))),
      accept = ord$accept, runtime_min = ord$runtime_min, mcmc = ord$settings$mcmc,
      waic = waic(ord$ll_resp)$waic, elpd = waic(ord$ll_resp)$elpd,
      holdout_ll = holdout_ll(ord, hold)$total, hit = hit_rate(ord, hold),
      lmd_nr = lmd_nr(ord$ll),
      ll_trace_halves = c(first = mean(ord$ll[seq_len(length(ord$ll) %/% 2)]),
                          second = mean(ord$ll[-seq_len(length(ord$ll) %/% 2)])))
    if (!is.null(ord2)) R$hb$seed_check <- list(
      max_abs_z = max(abs(pm(ord$draws$phibar) - pm(ord2$draws$phibar)) / psd(ord$draws$phibar)),
      cut_max_abs_z = max(abs(pm(ord$draws$cut) - pm(ord2$draws$cut)) / psd(ord$draws$cut)))
    if (!is.null(ordl)) R$hb$link <- list(
      waic_cloglog = R$hb$waic, waic_logit = waic(ordl$ll_resp)$waic,
      elpd_cloglog = R$hb$elpd, elpd_logit = waic(ordl$ll_resp)$elpd,
      holdout_cloglog = R$hb$holdout_ll, holdout_logit = holdout_ll(ordl, hold)$total,
      hit_cloglog = R$hb$hit, hit_logit = hit_rate(ordl, hold),
      lmd_nr_cloglog = R$hb$lmd_nr, lmd_nr_logit = lmd_nr(ordl$ll),
      phibar_logit = pm(ordl$draws$phibar), cut_logit = pm(ordl$draws$cut))
    if (!is.null(mnl)) R$hb$mnl <- list(phibar = pm(mnl$draws$phibar), sd = psd(mnl$draws$phibar),
                                        hit = hit_rate(mnl, hold))
  }

  # ---- D. hierarchical, respondent-specific cut points --------------------
  het <- load_fit("hetcut_cloglog", seg)
  if (!is.null(het) && !is.null(ord)) {
    P <- het$settings$P; W <- het$settings$W
    Sig <- apply(het$draws$Sigma, c(2, 3), mean)
    corr <- cov2cor(Sig)
    Phi_pm <- apply(het$draws$Phi, c(2, 3), mean)
    cuts_i <- phi_to_cutmat(Phi_pm, P, W)
    R$hetcut <- list(
      phibar = pm(het$draws$phibar)[1:P], phibar_sd = psd(het$draws$phibar)[1:P],
      cut_pop = pm(het$draws$phibar)[(P + 1):(P + W - 1)],
      sd_c1_across_resp = sqrt(Sig[P + 1, P + 1]),
      sd_logdiff_across_resp = sqrt(diag(Sig)[(P + 2):(P + W - 1)]),
      corr_beta_c1 = setNames(corr[1:P, P + 1], cn),
      cuts_i_quantiles = apply(cuts_i, 2, quantile, c(.05, .25, .5, .75, .95)),
      alpha_i_quantiles = apply(Fclog(cuts_i), 2, quantile, c(.05, .25, .5, .75, .95)),
      waic = waic(het$ll_resp)$waic, elpd = waic(het$ll_resp)$elpd,
      holdout_ll = holdout_ll(het, hold)$total, hit = hit_rate(het, hold),
      accept = het$accept, runtime_min = het$runtime_min)
  }

  # ---- E. practitioner move 1: dichotomise ---------------------------------
  t2bl <- load_fit("t2b_logit", seg); t2bc <- load_fit("t2b_cloglog", seg)
  topl <- load_fit("top_logit", seg); topc <- load_fit("top_cloglog", seg)
  k2l  <- load_fit("k2_logit", seg);  k3l  <- load_fit("k3_logit", seg)
  if (!is.null(ord) && !is.null(t2bl)) {
    sd_ratio <- function(f) psd(f$draws$phibar) / psd(ord$draws$phibar)
    bins <- list(t2b_logit = t2bl, t2b_cloglog = t2bc, top_logit = topl, top_cloglog = topc)
    bins <- bins[!vapply(bins, is.null, logical(1))]
    R$dichot <- list(
      sd_ratio = sapply(bins, sd_ratio),
      sd_ratio_mean = sapply(bins, function(f) mean(sd_ratio(f))),
      hit = c(ordinal = R$hb$hit, sapply(bins, hit_rate, dat_holdout = hold)),
      phibar = sapply(bins, function(f) pm(f$draws$phibar)))
    # the report's stepwise thresholds: four separate binary logits vs one
    # ordinal logit fit's four cut points
    steps <- list(k2 = k2l, k3 = k3l, k4 = t2bl, k5 = topl)
    if (all(!vapply(steps, is.null, logical(1))) && !is.null(ordl)) {
      R$dichot$stepwise <- data.frame(
        cut = paste0("c", 1:4),
        binary_mean = sapply(steps, function(f) mean(f$draws$cut[, 1])),
        binary_sd = sapply(steps, function(f) sd(f$draws$cut[, 1])),
        ordinal_mean = pm(ordl$draws$cut), ordinal_sd = psd(ordl$draws$cut),
        row.names = NULL)
      R$dichot$stepwise_beta_spread <- apply(sapply(steps, function(f) pm(f$draws$phibar)), 1,
                                             function(b) diff(range(b)))
    }
  }

  # ---- F. practitioner move 2: fixed probabilities -------------------------
  if (!is.null(ord)) {
    lw <- label_worth_tasks(ord, est)
    R$fixedprob <- list(
      incidence = heur_incidence(est$y),
      weights = wts_juster,
      label_worth = data.frame(label = 1:5,
        min = apply(lw, 2, min), q10 = apply(lw, 2, quantile, .10),
        median = apply(lw, 2, median),
        q90 = apply(lw, 2, quantile, .90), max = apply(lw, 2, max),
        fixed = wts_juster, row.names = NULL),
      label_worth_by_label_obs = sapply(1:5, function(w)
        quantile(lw[est$y == w, w], c(.1, .5, .9))))
  }

  # ---- G. managerial scenarios ---------------------------------------------
  if (!is.null(ord)) {
    Xs <- scenario_X(seg, cn)
    tab <- scenario_table(seg)
    hb_su <- function(f, k) if (is.null(f)) rep(NA, nrow(Xs)) else 100 * signup_hb(f, Xs, k)$mean
    hb_sd <- function(f, k) if (is.null(f)) rep(NA, nrow(Xs)) else 100 * signup_hb(f, Xs, k)$sd
    inc <- R$fixedprob$incidence
    sc <- data.frame(
      scenario = tab$scenario,
      report = tab$report_signup,
      ordinal_t2b = hb_su(ord, 4), ordinal_t2b_sd = hb_sd(ord, 4),
      ordinal_logit_t2b = hb_su(ordl, 4),
      binary_t2b_logit = hb_su(t2bl, 2), binary_t2b_logit_sd = hb_sd(t2bl, 2),
      binary_t2b_cloglog = hb_su(t2bc, 2),
      ordinal_top = hb_su(ord, 5), binary_top_logit = hb_su(topl, 2),
      expected_p = 100 * expected_p_hb(ord, Xs)$mean,
      heur_juster = 100 * inc["juster"], heur_t2b = 100 * inc["t2b"],
      row.names = NULL)
    if (!is.null(het)) sc$hetcut_t2b <- hb_su(het, 4)
    R$scenarios <- list(table = sc, levels = tab[, c("scenario", paste0("att", 1:9))])
  }

  # ---- H. specification tests at the hierarchical level --------------------
  exl <- load_fit("ext_lambda", seg); exd <- load_fit("ext_delta", seg)
  th <- function(f, null) if (is.null(f)) NULL else {
    t <- f$draws$theta
    list(mean = mean(t), sd = sd(t), ci = quantile(t, c(.025, .975)), null = null,
         prob_below_null = mean(t < null), accept = f$accept$theta,
         z = (mean(t) - null) / sd(t)) }
  exlh <- load_fit("ext_lambda_hetcut", seg); exdh <- load_fit("ext_delta_hetcut", seg)
  R$hb_spec <- list(lambda = th(exl, 1), delta = th(exd, 0),
                    lambda_hetcut = th(exlh, 1), delta_hetcut = th(exdh, 0))

  # Parametric bootstrap of the lambda tests: data simulated from the fitted
  # hierarchical model with lambda = 1 exactly (individual tastes drawn from
  # the estimated population distribution, so the simulated panel carries the
  # fitted heterogeneity), then both the aggregate and the hierarchical lambda
  # estimates. The real-data estimates are read against these. The expensive
  # hierarchical refits come from a separate script and are cached; the
  # aggregate part is recomputed here if that cache is absent.
  if (!is.null(ord)) {
    bf <- file.path(FD, sprintf("boot_lambda_%s.rds", seg))
    if (file.exists(bf)) {
      bt <- readRDS(bf)
      R$attenuation <- list(boot = bt,
        lambda_agg_sim_mean = mean(bt$lambda_agg), lambda_hb_sim_mean = mean(bt$lambda_hb))
    } else {
      phibar_pm <- pm(ord$draws$phibar); Sig_pm <- apply(ord$draws$Sigma, c(2, 3), mean)
      cutv <- pm(ord$draws$cut)
      lam_sim <- vapply(1:5, function(s) {
        B <- draw_betas(est$N, phibar_pm, Sig_pm, seed = 1000 + s)
        sim <- simulate_dual_hb(est$design, B, cutv, link = "cloglog", seed = 2000 + s)
        m <- fit_dual_mle(sim)
        fit_dual_mle_ext(sim, "lambda", m$par, 1)$theta
      }, numeric(1))
      R$attenuation <- list(lambda_agg_sim = lam_sim, lambda_agg_sim_mean = mean(lam_sim))
    }
    R$attenuation$lambda_agg_real <- a$spec$estimate[a$spec$test == "lambda = 1"]
    R$attenuation$lambda_hb_real <- if (!is.null(exl)) mean(exl$draws$theta) else NA
    bh <- file.path(FD, sprintf("boot_lambda_hetcut_%s.rds", seg))
    if (file.exists(bh)) {
      bt2 <- readRDS(bh)
      R$attenuation$boot_hetcut <- bt2
      R$attenuation$lambda_hetcut_sim_mean <- mean(bt2$lambda_hb)
    }
    R$attenuation$lambda_hetcut_real <- if (!is.null(exlh)) mean(exlh$draws$theta) else NA
  }

  results[[seg]] <- R
}

results$settings <- list(mcmc = list(R = 30000, burn = 12000, thin = 10),
                         n_holdout = 2L, W = 5L, J = 3L, T = 12L,
                         wts_juster = wts_juster, date = as.character(Sys.Date()))
saveRDS(results, file.path(OUT, "empirical_cobrand.rds"))
message("wrote ", file.path(OUT, "empirical_cobrand.rds"))

# ---- figures -----------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
seg_lab <- c(member = "Members", nonmember = "Non-members")

# label worth: the purchase probability behind each scale point across tasks
lwd <- do.call(rbind, lapply(c("member", "nonmember"), function(s) {
  d <- results[[s]]$fixedprob$label_worth; if (is.null(d)) return(NULL)
  d$segment <- seg_lab[s]; d }))
if (!is.null(lwd)) {
  p <- ggplot(lwd, aes(x = factor(label))) +
    geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 2.2, colour = "#317EA7", alpha = .5) +
    geom_linerange(aes(ymin = min, ymax = max), linewidth = .6, colour = "#317EA7") +
    geom_point(aes(y = median), colour = "#317EA7", size = 2) +
    geom_point(aes(y = fixed), colour = "#BA403D", shape = 4, size = 3, stroke = 1.2) +
    facet_wrap(~ segment) +
    labs(x = "Scale point", y = "Purchase probability behind the label",
         title = NULL) +
    scale_y_continuous(limits = c(0, 1)) + theme_minimal(base_size = 11)
  ggsave(file.path(FIG, "empirical_label_worth.png"), p, width = 6.5, height = 3.2, dpi = 300)
}

# scenarios: sign-up by scenario and method
scd <- do.call(rbind, lapply(c("member", "nonmember"), function(s) {
  d <- results[[s]]$scenarios$table; if (is.null(d)) return(NULL)
  long <- rbind(
    data.frame(scenario = d$scenario, method = "Ordinal model (top-two-box)", value = d$ordinal_t2b),
    data.frame(scenario = d$scenario, method = "Binary top-two-box (report's model)", value = d$binary_t2b_logit),
    data.frame(scenario = d$scenario, method = "Fixed probabilities", value = d$heur_juster),
    data.frame(scenario = d$scenario, method = "Report", value = d$report))
  long$segment <- seg_lab[s]; long }))
if (!is.null(scd)) {
  scd$scenario <- factor(scd$scenario, levels = c("Current", paste0("S", 1:7), "Worst"))
  p <- ggplot(scd, aes(scenario, value, colour = method, group = method)) +
    geom_line() + geom_point(size = 1.8) + facet_wrap(~ segment) +
    labs(x = NULL, y = "Likely to sign up (%)", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(FIG, "empirical_scenarios.png"), p, width = 7, height = 3.6, dpi = 300)
}
message("figures written to ", FIG)
