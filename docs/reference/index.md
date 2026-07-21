# Package index

## Main workflow

Functions most users need for an analysis.

- [`dr_gpsm()`](https://leo-liuqiang.github.io/BC-GPSM/reference/dr_gpsm.md)
  : Doubly-robust GPS matching estimator
- [`gps_pre_process()`](https://leo-liuqiang.github.io/BC-GPSM/reference/gps_pre_process.md)
  : Generalized Propensity Score (GPS) Pre-processing (with optional GPS
  fitting)
- [`gps_histogram()`](https://leo-liuqiang.github.io/BC-GPSM/reference/gps_histogram.md)
  : Plot overlaid histograms of non-crossfitted diagnostic GPS
- [`balance_check_plot()`](https://leo-liuqiang.github.io/BC-GPSM/reference/balance_check_plot.md)
  : Covariate balance check plot

## Matching and contrasts

Lower-level tools used by the main estimator.

- [`gps_matching()`](https://leo-liuqiang.github.io/BC-GPSM/reference/gps_matching.md)
  : GPS- or Covariate-matching ATE estimator (numeric-column interface)
- [`build_contrast()`](https://leo-liuqiang.github.io/BC-GPSM/reference/build_contrast.md)
  : Create all-pairwise contrast matrix

## Propensity-score helpers

Advanced helpers for cross-fitting and GPS log-ratio construction.

- [`pmatch_fit_ps_fold()`](https://leo-liuqiang.github.io/BC-GPSM/reference/pmatch_fit_ps_fold.md)
  : Fit PS model on training fold and predict probabilities on eval fold
- [`pmatch_compute_V()`](https://leo-liuqiang.github.io/BC-GPSM/reference/pmatch_compute_V.md)
  : Compute log-ratio PS index V from multinomial probabilities
