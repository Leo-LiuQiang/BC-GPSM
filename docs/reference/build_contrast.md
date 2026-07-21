# Create all-pairwise contrast matrix

Create all-pairwise contrast matrix

## Usage

``` r
build_contrast(levels, ref = NULL)
```

## Arguments

- levels:

  Character vector of treatment levels

- ref:

  Reference treatment level. If `NULL`, the last level is used.

## Value

A contrast matrix of size Choose(K, 2) x K with entries -1, 0, or 1

## Examples

``` r
build_contrast(c("A", "B", "C"), ref = "A")
#>      A  B C
#> BvA -1  1 0
#> CvA -1  0 1
#> CvB  0 -1 1
build_contrast(c("A", "B", "C"))
#>      A B  C
#> BvA -1 1  0
#> AvC  1 0 -1
#> BvC  0 1 -1
```
