# BC.GPSM 0.1.0

## Initial CRAN release

* Implements bias-corrected generalized propensity score matching for binary
  and multi-valued treatments.
* Supports explicit or default reference treatments, nearest-neighbor matching
  with configurable matching ratios, outcome adjustment, two-step calibration,
  and bootstrap confidence intervals.
* Provides generalized propensity-score overlap and covariate-balance
  diagnostics, including publication-ready love plots.
* Provides multinomial logistic GPS estimation and linear outcome adjustment as
  the baseline workflow, with optional GAM, GBM, random forest, ranger, and
  XGBoost learners.
* Treats the ranger and XGBoost integrations as experimental optional backends
  while their statistical behavior is evaluated further.
* Includes a runnable package vignette, pkgdown documentation, actionable input
  diagnostics, and reliability tests across macOS, Windows, and Linux.
