# Fit PS model on training fold and predict probabilities on eval fold

Fit PS model on training fold and predict probabilities on eval fold

## Usage

``` r
pmatch_fit_ps_fold(
  train_df,
  eval_df,
  t_name,
  x_names,
  model = c("logit", "gbm", "gam", "xgboost", "ranger"),
  tune = FALSE,
  tune_control = NULL,
  tune_grid = NULL,
  gps_params = NULL,
  tr_levels = NULL,
  seed = 12345
)
```

## Arguments

- train_df:

  training data frame

- eval_df:

  evaluation data frame

- t_name:

  treatment column name (factor)

- x_names:

  covariate column names

- model:

  "logit","gbm","gam","xgboost","ranger". Flexible learners require
  optional packages; XGBoost and ranger are experimental.

- tune:

  Logical; optional caret tuning for logit/gbm.

- tune_control:

  [`caret::trainControl()`](https://rdrr.io/pkg/caret/man/trainControl.html)
  or NULL.

- tune_grid:

  named list with components logit/gbm (same style as your
  gps_tune_grids)

- gps_params:

  optional named list of model params (same style as your gps_params)

- tr_levels:

  full set of treatment levels to enforce (character)

- seed:

  integer seed

## Value

list(e_hat_eval = n_eval x K matrix, levels = tr_levels)

## Examples

``` r
set.seed(1)
n <- 75
dat <- data.frame(
  trt = factor(rep(c("A", "B", "C"), each = 25)),
  x1 = rnorm(n),
  x2 = runif(n)
)

train_id <- c(1:20, 26:45, 51:70)
ps_fit <- pmatch_fit_ps_fold(
  train_df = dat[train_id, ],
  eval_df = dat[-train_id, ],
  t_name = "trt",
  x_names = c("x1", "x2"),
  model = "logit",
  tr_levels = levels(dat$trt)
)

head(ps_fit$e_hat_eval)
#>            A         B         C
#> 21 0.3063262 0.3083768 0.3852970
#> 22 0.2382435 0.3962101 0.3655464
#> 23 0.2211287 0.4612795 0.3175918
#> 24 0.3320404 0.4342977 0.2336620
#> 25 0.3891431 0.2445209 0.3663360
#> 46 0.2911563 0.4175741 0.2912696
```
