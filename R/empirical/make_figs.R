# Two figures for Section 4 that need only the blind summary, so they can be
# redrawn without the client data:
#   empirical_thresholds.png  cut points on the probability scale, alpha_w = F(c_w),
#                             aggregate ML (clustered SE) and common-cut HB
#   empirical_cuts_i.png      the spread of respondent-specific cut points across
#                             people, utility scale, with the common-cut estimate
# The label-worth and scenario figures are written by analyze_cobrand.R.
#   Rscript R/empirical/make_figs.R      (from the repo root)
suppressPackageStartupMessages(library(ggplot2))
res <- readRDS(file.path("R", "output", "empirical_cobrand.rds"))
FIG <- file.path("R", "output", "figs")
seg_lab <- c(member = "Members", nonmember = "Non-members")
boundaries <- c("very unlikely | somewhat unlikely", "somewhat unlikely | neither",
                "neither | somewhat likely", "somewhat likely | very likely")
col_est <- c("Aggregate ML" = "#616D64", "HB, common cut points" = "#317EA7",
             "HB, respondent-specific cut points" = "#B28000")
th <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

# --- thresholds on the probability scale -------------------------------------
a <- do.call(rbind, lapply(names(seg_lab), function(s) { R <- res[[s]]
  rbind(data.frame(segment = seg_lab[[s]], est = names(col_est)[1], w = 1:4, alpha = R$agg$alpha,
                   lo = pmax(0, R$agg$alpha - 1.96 * R$agg$se_alpha_cl),
                   hi = pmin(1, R$agg$alpha + 1.96 * R$agg$se_alpha_cl)),
        data.frame(segment = seg_lab[[s]], est = names(col_est)[2], w = 1:4, alpha = R$hb$alpha,
                   lo = R$hb$alpha_ci[1, ], hi = R$hb$alpha_ci[2, ])) }))
a$boundary <- factor(a$w, labels = boundaries)
a$est <- factor(a$est, levels = names(col_est)[1:2])
p <- ggplot(a, aes(x = alpha, y = boundary, colour = est)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .5),
                  size = .4, linewidth = .6) +
  facet_wrap(~ segment) +
  scale_colour_manual(values = col_est[1:2], name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(x = "Purchase probability at which the label begins", y = NULL) + th
ggsave(file.path(FIG, "empirical_thresholds.png"), p, width = 7, height = 3, dpi = 300)

# --- respondent-specific cut points: the spread across people ----------------
# Utility scale: on the probability scale most respondents' alpha_i sit at 0
# or 1 and the picture is two walls.
ci <- do.call(rbind, lapply(names(seg_lab), function(s) { R <- res[[s]]
  q <- R$hetcut$cuts_i_quantiles
  data.frame(segment = seg_lab[[s]], w = 1:4, q05 = q["5%", ], q25 = q["25%", ], q50 = q["50%", ],
             q75 = q["75%", ], q95 = q["95%", ], common = R$hb$cut) }))
ci$boundary <- factor(ci$w, labels = boundaries)
lab_i <- "respondent-specific: median, 25-75% and 5-95% of respondents"
lab_c <- "common cut point"
p2 <- ggplot(ci, aes(y = boundary)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .4) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = .6, colour = col_est[[3]]) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2.2, colour = col_est[[3]], alpha = .5) +
  geom_point(aes(x = q50, colour = lab_i), size = 2) +
  geom_point(aes(x = common, colour = lab_c), shape = 18, size = 3.2) +
  facet_wrap(~ segment) +
  scale_colour_manual(values = setNames(c(col_est[[2]], col_est[[3]]), c(lab_c, lab_i)), name = NULL,
                      guide = guide_legend(nrow = 2)) +
  labs(x = "Cut point on the utility scale", y = NULL) + th
ggsave(file.path(FIG, "empirical_cuts_i.png"), p2, width = 7, height = 3.2, dpi = 300)
cat("wrote", file.path(FIG, c("empirical_thresholds.png", "empirical_cuts_i.png")), sep = "\n")
