
## Persistent outside option (opened 2026-10-06)

Status: notebook written and smoke-tested; **not yet rendered at full scale**.

- [ ] Render `notebooks/09-persistent-outside-option.qmd` (~16-18 min, 48 HB fits
      at N=300, T=12, ITER=15000). Writes `R/output/persistent_outside.rds`.
- [ ] Read the rendered numbers back into @sec-persist in `sections/02_model.qmd`
      — the prose currently asserts the mechanism but cites no estimate.
- [ ] Decide whether the sampler variant belongs in `R/hierarchical/hb.R` as a
      `persist` flag, or stays notebook-local. Notebook 09 inlines a trimmed
      sampler (no het_cut, no choice_only) rather than extending the library one.
- [ ] Apply the sigma_z^2 vs pi^2/6 test to the empirical data in Section 4.

Verified at full scale before writing (single fit, N=300, T=12, ITER=15000):
  sigma_z^2 posterior mean 1.623, 95% CI [1.254, 2.061], target pi^2/6 = 1.645.
  beta-bar RMSE 0.0499. One fit = 20.2 sec.

Known wrinkle, already handled in the notebook prose: at a true sigma_z^2 = 0 the
estimate returns ~0.22, not 0. Each z_i is estimated from T answers so carries
noise, and the inverse-gamma posterior has no density at the boundary. That row
is the design's detection floor, not a coverage failure — do not let it get
re-reported as "the test is biased."
