# Outside Good Uncertainty

## Ordinal Dual Response in Choice-Based Conjoint Analysis

In choice-based conjoint studies, respondents often first pick their preferred product and then state whether they would actually buy it, a design known as the dual-response format. We provide a behavioral foundation for this format and extend the follow-up purchase question from a yes/no answer to an ordinal scale (e.g., "definitely not" to "definitely would buy"). In our model, consumers identify their most-preferred option among a set of alternatives but face uncertainty in their evaluation of the outside good. We develop estimation routines for both aggregate and hierarchical Bayesian versions of the model and apply them to a conjoint study conducted with a consumer insights consultancy. The model outperforms two common alternatives: collapsing the ordinal scale into a binary response, and assigning a fixed purchase probability to each scale point.


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
and must never be committed. It also holds the ingest script that reads that
data and every cached fit (`derived/`); the tracked side of the application is
`R/empirical/`, and the only artefact that crosses from one to the other is the
blind summary `R/output/empirical_cobrand.rds`. See `R/README.md`. Anything under
`conference/` other than `sawtooth-2027/` is local material and is ignored by
default.
