
## Persistent outside option — CLOSED 2026-10-08

Cut from the paper. The pi^2/6 prediction depended on the respondent having
resolved the outside-good shock, and that story was ruled out for ordinal
follow-ups (a resolved respondent has a 0/1 purchase probability and only the
endpoints of a likelihood scale to report). Without it the extension reduces to
unrestricted location heterogeneity, which sec-hetcut already covers.

Notebook moved to archive/09-persistent-outside-option.qmd; results file and
freeze removed. Recoverable from commit 083bb09 if the framing is ever revived.

## Missing PDFs for cited work — OPEN 2026-11-02

Seven works the paper cites have no PDF in `lit/`. Download and file them:

- Sawtooth Software (2017), *The CBC System for Choice-Based Conjoint Analysis, Version 9*, technical paper (sawtoothsoftware.com) → `lit/1 dual response and the no-choice option/`
- King, Murray, Salomon & Tandon (2004), anchoring vignettes, *APSR* 98(1) → `lit/2 discrete choice models and methods/`
- Yatchew & Griliches (1985), Specification Error in Probit Models, *REStat* 67(1) → `lit/2 discrete choice models and methods/`
- Ben-Akiva & Morikawa (1990), switching models from RP and SI data, *Transportation* 17 → `lit/2 discrete choice models and methods/`
- Swait & Louviere (1993), The Role of the Scale Parameter, *JMR* 30(3) → `lit/2 discrete choice models and methods/`
- de Palma & Kilani (2007), Invariance of Conditional Maximum Utility, *JET* 132(1), 137–146 (post-print on HAL, hal-04826236) → `lit/2 discrete choice models and methods/`; cited in Appendices A and D beside McFadden (1981)
- Agresti (2002), *Categorical Data Analysis*, 2nd ed., Wiley (book; a chapter scan of the ordinal-regression chapter is enough) → `lit/2 discrete choice models and methods/`; cited with the Christensen ordinal vignette (filed 2026-11-03) for the cumulative-link name and the link naming footnote in §2

Note: the Sawtooth CBC/HB paper filed on 2026-11-02 is V5.6 (2021); the bib
cites the v5 paper of April 2016. Same series, later revision; check the
dual-response coding passage still matches before relying on it.

## The outside-good normalization and the probability-scale quantities — OPEN 2026-11-03

$V_{it0} = 0$ is a normalization for identifying $(\beta, \mathbf{c})$ (Section 2
says so), but the probability-scale quantities are not invariant to it.
Re-coding the reference levels shifts every $V(x)$ by a constant $k$; the fitted
cut points shift by $k$ to keep the likelihood, so $\alpha_w = F(c_w)$,
$E[p] = S/(S+1)$ and the label-worth values $E[p \mid y = w, S]$ all change.
The data identify $c_w + k$, not $c_w$: the outside good's location and the
thresholds' level trade off, and the model resolves it by placing the outside
good at the utility of the all-reference-levels card. To do:

- State this in the identification section (one paragraph after "The
  outside good's location"): which quantities are coding-invariant (label
  shares $\Pr(y \ge w \mid x)$, model comparisons at a margin, the $\lambda$ and
  $\zeta$ tests, threshold spacings, the set-dependence of label worth) and
  which are relative to the normalization (individual $\alpha_w$, $E[p]$,
  individual label-worth values).
- Reword Section 4 where $\alpha_w$, $E[p]$ and label worth are read as
  absolute probabilities ("somewhat likely begins at 0.57", the 73% / 82%
  take rates, "the constants miss every label"): keep the spacing and
  set-dependence claims, qualify the levels.
- Add to the limitations list: pinning the level needs a calibrated outside
  good (observed sign-ups) or a calibrated label; the binary AnyMax model has
  the mirror-image assumption (label fixed, constant estimated).
- Consider: with sign-up data for even one card, $k$ is identified and the
  probability scale becomes absolute; note as the natural validation study.

Also still to do from the same session: Section 4 should adopt the
calibrated $\lambda$ (0.46 [0.30, 0.62] members, 0.47 [0.27, 0.68]
non-members, `calib` fields in the blind summary) in place of the
"net of the bootstrap" narrative; and the dichotomisation paragraph should
say that the 22% / 18% against the incumbent logit mixes link and
dichotomisation (same-link cost 9% / 8% at top-two-box, 14% / 19% at top-box,
`dichot$sd_ratio_mean` and the new `dichot$sd_ratio_ind_mean`).
