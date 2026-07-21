# Doubly-robust GPS matching estimator

Implements a doubly-robust generalized propensity score (GPS) matching
estimator with optional cross-fitting and bootstrap inference. GPS can
be estimated via multinomial logistic regression, GBM, GAM, XGBoost, or
ranger; outcome models can be fitted via LM, RF, GBM, GAM, XGBoost, or
ranger. Flexible learners and tuning integrations require optional
packages. The XGBoost and ranger backends are experimental and remain
under active testing.

## Usage

``` r
dr_gpsm(
  data,
  treatment,
  treatment_ref = NULL,
  covariate,
  outcome,
  gps_model = c("logit", "gbm", "gam", "xgboost", "ranger"),
  outcome_model = c("none", "lm", "rf", "gbm", "gam", "xgboost", "ranger"),
  folds = 2,
  fold_seed = 12345,
  nboot = 500,
  boot_weight = c("multinom", "exp"),
  hist = FALSE,
  hist_path = NULL,
  gps_params = NULL,
  outcome_params = NULL,
  match_on = c("gps", "covariates"),
  cov_distance = c("mahalanobis", "euclidean"),
  standardize = FALSE,
  match_ratio = 1L,
  gps_tune = FALSE,
  gps_tune_control = NULL,
  gps_tune_grids = NULL,
  gps_seed = 12345,
  outcome_tune = FALSE,
  outcome_tune_control = NULL,
  outcome_tune_grids = NULL,
  outcome_seed = 12345,
  two_step_calibration = FALSE,
  calib_shrinkage = 1
)
```

## Arguments

- data:

  Data frame with treatment, covariates, and outcome variables

- treatment:

  Column index of the treatment variable

- treatment_ref:

  Reference treatment level. If `NULL`, the last observed treatment
  level is used. The selected reference is moved to the first factor
  level internally for GPS log-ratio construction.

- covariate:

  Numeric vector of covariate column indices

- outcome:

  Column index of the outcome variable

- gps_model:

  Choice of GPS model: 'logit', 'gbm', 'gam', 'xgboost', or 'ranger'.
  All choices except 'logit' require optional packages.

- outcome_model:

  Choice of outcome model: 'none','lm','rf','gbm','gam', 'xgboost', or
  'ranger'. Flexible learners require optional packages.

- folds:

  Number of folds for cross-fitting (default = 2)

- fold_seed:

  Integer seed for fold assignment (default 12345).

- nboot:

  Number of bootstrap replications (default = 500)

- boot_weight:

  Bootstrap reweighting scheme: 'multinom' or 'exp'.

- hist:

  Logical; TRUE to save a diagnostic full-sample GPS histogram

- hist_path:

  File path (pdf/png/jpg) to save histogram if hist = TRUE

- gps_params:

  Optional named list of hyperparameters per GPS model

- outcome_params:

  Optional named list of hyperparameters per outcome model

- match_on:

  'gps' (default) or 'covariates'

- cov_distance:

  For covariate matching: 'euclidean' or 'mahalanobis' (default
  'mahalanobis')

- standardize:

  If TRUE, center/scale covariate features before matching (default
  FALSE)

- match_ratio:

  Integer \>= 1; number of matches per unit per target group (default 1)

- gps_tune:

  Logical; if TRUE, enable optional caret tuning for GPS
  ('logit','gbm').

- gps_tune_control:

  Optional
  [`caret::trainControl()`](https://rdrr.io/pkg/caret/man/trainControl.html)
  for GPS tuning; default is 5-fold CV

- gps_tune_grids:

  Optional named list of GPS tuning grids (e.g.,
  `list(gbm = ..., logit = ...)`)

- gps_seed:

  Integer seed for GPS tuning and fitting (default 12345)

- outcome_tune:

  Logical; if TRUE, enable optional caret tuning for outcome models
  ('rf','gbm').

- outcome_tune_control:

  Optional
  [`caret::trainControl()`](https://rdrr.io/pkg/caret/man/trainControl.html)
  for outcome tuning

- outcome_tune_grids:

  Optional named list of outcome tuning grids (e.g.,
  `list(rf = ..., gbm = ...)`)

- outcome_seed:

  Integer seed for outcome tuning and fitting (default 12345)

- two_step_calibration:

  Logical; if TRUE, apply two-step residual calibration using ridge on
  GPS index V.

- calib_shrinkage:

  Numeric in `[0,1]`; shrinkage multiplier applied to the calibration
  step size gamma.

## Value

A list with components:

- estimate:

  Named numeric vector of doubly-robust ATE estimates for each contrast.

- ci_lower:

  Named numeric vector of lower 95% confidence bounds.

- ci_upper:

  Named numeric vector of upper 95% confidence bounds.

## Examples

``` r
set.seed(1)
n <- 75
dat <- data.frame(
  trt = factor(rep(c("A", "B", "C"), each = 25)),
  x1 = rnorm(n),
  x2 = runif(n)
)
dat$y <- 1 + 0.4 * (dat$trt == "B") + 0.8 * (dat$trt == "C") +
  dat$x1 + rnorm(n, sd = 0.5)

fit <- dr_gpsm(
  data = dat,
  treatment = 1,
  treatment_ref = "A",
  covariate = 2:3,
  outcome = 4,
  gps_model = "logit",
  outcome_model = "lm",
  folds = 2,
  nboot = 10
)

fit$estimate
#>       BvA       CvA       CvB 
#> 0.5887164 0.9360901 0.3473737 
```
