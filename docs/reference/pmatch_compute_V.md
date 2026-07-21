# Compute log-ratio PS index V from multinomial probabilities

Compute log-ratio PS index V from multinomial probabilities

## Usage

``` r
pmatch_compute_V(e_mat, ref_col = ncol(e_mat), eps = 1e-12)
```

## Arguments

- e_mat:

  n x K matrix of probabilities (rows sum to 1)

- ref_col:

  Column index or column name of the reference treatment. Defaults to
  the last column.

- eps:

  small clipping value

## Value

n x (K-1) matrix of log-ratios relative to the reference column.

## Examples

``` r
e_mat <- matrix(
  c(0.20, 0.50, 0.30,
    0.40, 0.40, 0.20),
  nrow = 2,
  byrow = TRUE
)
pmatch_compute_V(e_mat)
#>            [,1]      [,2]
#> [1,] -0.4054651 0.5108256
#> [2,]  0.6931472 0.6931472
pmatch_compute_V(e_mat, ref_col = 1)
#>           [,1]       [,2]
#> [1,] 0.9162907  0.4054651
#> [2,] 0.0000000 -0.6931472
```
