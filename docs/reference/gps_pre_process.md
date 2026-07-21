# Generalized Propensity Score (GPS) Pre-processing (with optional GPS fitting)

In the cross-fitted BC-GPSM workflow, this function should primarily:

1.  enforce treatment factor + reference level;

2.  keep only complete cases (treatment + covariates);

3.  ALWAYS create gps\_ and loggps\_ placeholder columns (NA by
    default);

4.  optionally fit a full-sample GPS model (diagnostic only, e.g., for
    histogram).

## Usage

``` r
gps_pre_process(
  data,
  treatment,
  treatment_ref = NULL,
  covariate,
  gps_model = c("logit", "gbm", "gam", "xgboost", "ranger"),
  gps_params = NULL,
  tune = FALSE,
  tune_control = NULL,
  tune_grids = NULL,
  seed = 12345,
  fit_gps = TRUE
)
```

## Arguments

- data:

  A data.frame containing treatment/covariates (and possibly outcome).

- treatment:

  Single numeric column index for treatment.

- treatment_ref:

  Optional reference level. If `NULL`, the last observed treatment level
  is used. The selected reference is moved to the first factor level
  internally because the log-GPS ratios use it as the denominator.

- covariate:

  Numeric vector of covariate column indices.

- gps_model:

  One of "logit","gbm","gam","xgboost","ranger" (used only if
  fit_gps=TRUE). Flexible learners require optional packages; the
  XGBoost and ranger backends are experimental and under active testing.

- gps_params:

  Optional named list of model-specific parameters.

- tune:

  Logical; if TRUE, optional caret tuning is used for "logit" and "gbm"
  when fit_gps=TRUE.

- tune_control:

  Optional
  [`caret::trainControl()`](https://rdrr.io/pkg/caret/man/trainControl.html)
  object for tuning.

- tune_grids:

  Optional tuning grids list.

- seed:

  Integer seed.

- fit_gps:

  Logical; if FALSE, do NOT fit GPS; return placeholders only.

## Value

data.frame with original columns + gps\_ and loggps\_ (always present).

## Examples

``` r
set.seed(1)
n <- 75
dat <- data.frame(
  trt = factor(rep(c("A", "B", "C"), each = 25)),
  x1 = rnorm(n),
  x2 = runif(n)
)

gps_dat <- gps_pre_process(
  data = dat,
  treatment = 1,
  covariate = 2:3,
  gps_model = "logit"
)

levels(gps_dat$trt)
#> [1] "C" "A" "B"
head(gps_dat[, grep("^(gps|loggps)_", names(gps_dat))])
#>       gps_C     gps_A     gps_B     loggps_A     loggps_B
#> 1 0.3144419 0.2994923 0.3860658 -0.048710823  0.205208285
#> 2 0.3361770 0.3288469 0.3349761 -0.022045534 -0.003578721
#> 3 0.3381825 0.3544560 0.3073615  0.046998320 -0.095561236
#> 4 0.3682619 0.3778625 0.2538756  0.025735941 -0.371950044
#> 5 0.3439777 0.3434879 0.3125344 -0.001424788 -0.095862102
#> 6 0.3491780 0.3849738 0.2658482  0.097593494 -0.272656556
```
