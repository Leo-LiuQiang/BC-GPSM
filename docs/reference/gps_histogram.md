# Plot overlaid histograms of non-crossfitted diagnostic GPS

Plot overlaid histograms of non-crossfitted diagnostic GPS

## Usage

``` r
gps_histogram(data, bins = 30, palette = NULL, alpha = 0.4, eps = 1e-06)
```

## Arguments

- data:

  Data frame returned by
  [`gps_pre_process()`](https://leo-liuqiang.github.io/BC-GPSM/reference/gps_pre_process.md)
  (must contain `gps_` columns)

- bins:

  Number of histogram bins (default = 30)

- palette:

  Character vector of fill colors (recycled if shorter than the number
  of GPS columns)

- alpha:

  Fill transparency (0-1, default = 0.4) for overlapping areas

- eps:

  Small value to truncate probabilities (avoids Inf after logit; default
  = 1e-6)

## Value

A ggplot object showing overlaid histograms and density curves of GPS
distributions

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
  treatment_ref = "A",
  covariate = 2:3,
  gps_model = "logit"
)

gps_histogram(gps_dat, bins = 10)
```
