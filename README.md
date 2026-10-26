# Outside Good Uncertainty

## Ordinal Dual Response in Choice-Based Conjoint Analysis

We present a behavioral model in which consumers identify their most-preferred option among a set of alternatives but face uncertainty in their evaluation of the outside good. The model rationalizes the dual-response format used in choice-based conjoint analysis and extends the follow-up purchase question from binary to ordinal. The resulting likelihood function combines a multinomial logit for the forced choice with a cumulative-link ordinal model driven by the task's inclusive value. The link is logistic or Gumbel, depending on whether the respondent resolves her uncertainty before answering, and the binary dual-response model is a special case. We derive how much information a finer response scale adds, provide estimation routines for both aggregate and hierarchical Bayesian versions of the model, and show that the model outperforms both dichotomizing the scale into a binary response and assigning fixed purchase probabilities to each scale point. We apply the model to a choice-based conjoint study conducted in partnership with a consumer insights consultancy.


## Repository layout

`lit/` is a **symlink into Dropbox**, not tracked content:

| path in this project | real location |
|---|---|
| `lit/` | `~/Dropbox/research_projects/2026_ordinal_dualresponse/lit` |

It holds publisher-copyright PDFs, so it stays out of git but lives in Dropbox to
remain backed up and shareable. A fresh clone will not have it; recreate it with
`./setup-links.sh`.

Note the `.gitignore` entry is written **without a trailing slash**. Git treats a
symlink as a file rather than a directory, so `lit/` would silently fail to match
once `lit` became a link. Do not "fix" it by putting the slash back.

`empirical application/` is likewise untracked: it holds un-anonymized client data
and must never be committed. Anything under `conference/` other than
`sawtooth-2027/` is local material and is ignored by default.
