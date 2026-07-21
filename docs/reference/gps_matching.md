# GPS- or Covariate-matching ATE estimator (numeric-column interface)

For cross-fitted inference, use dr_gpsm(); gps_matching()'s bootstrap is
intended for standalone (non-cross-fitted) usage.

## Usage

``` r
gps_matching(
  data,
  treatment,
  outcome,
  pred,
  contrast,
  nboot = NULL,
  match_on = c("gps", "covariates"),
  covariate = NULL,
  cov_distance = c("mahalanobis", "euclidean"),
  standardize = FALSE,
  ridge = 1e-08,
  match_ratio = 1L,
  return_tau = FALSE,
  do_boot = FALSE,
  boot_weight = c("multinom", "exp")
)
```

## Arguments

- data:

  Data frame

- treatment:

  Column index of the treatment variable (factor)

- outcome:

  Column index of the outcome variable

- pred:

  pred K x n matrix of m_k(X) predictions aligned with rows of data

- contrast:

  C(K, 2) x K contrast matrix from build_contrast()

- nboot:

  Number of bootstrap replications (only used if do_boot=TRUE; can be
  NULL otherwise)

- match_on:

  'gps' (default) or 'covariates'

- covariate:

  Numeric vector of covariate column indices (required if
  match_on='covariates')

- cov_distance:

  'euclidean' or 'mahalanobis' (only for match_on='covariates')

- standardize:

  Logical; standardize covariate features before matching

- ridge:

  Small ridge for covariance in Mahalanobis whitening

- match_ratio:

  Integer \>= 1; number of matches per unit per target group

- return_tau:

  Logical; if TRUE, also return tau (and tau_centered)

- do_boot:

  Logical; if TRUE, compute bootstrap CI within this fold; if FALSE,
  skip (ci\_\* = NA)

- boot_weight:

  Bootstrap reweighting scheme: 'multinom' or 'exp' (only used when
  `do_boot=TRUE`).

## Value

list(estimate, ci_lower, ci_upper, tau?, tau_centered?)

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

gps_dat <- gps_pre_process(
  dat,
  treatment = 1,
  treatment_ref = "A",
  covariate = 2:3,
  gps_model = "logit"
)
contrast <- build_contrast(levels(gps_dat$trt), ref = "A")
arm_means <- tapply(gps_dat$y, gps_dat$trt, mean)[levels(gps_dat$trt)]
pred <- matrix(rep(as.numeric(arm_means), each = nrow(gps_dat)),
  nrow = length(arm_means),
  byrow = TRUE
)

gps_matching(
  data = gps_dat,
  treatment = 1,
  outcome = 4,
  pred = pred,
  contrast = contrast,
  do_boot = FALSE
)
#> $estimate
#>       BvA       CvA       CvB 
#> 0.3063032 0.8779711 0.5716679 
#> 
#> $ci_lower
#> BvA CvA CvB 
#>  NA  NA  NA 
#> 
#> $ci_upper
#> BvA CvA CvB 
#>  NA  NA  NA 
#> 
```
