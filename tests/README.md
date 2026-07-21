# Running BC.GPSM Tests

Run commands from the package root.

## Focused Reliability Tests

```sh
Rscript -e 'devtools::test(filter = "dr-gpsm-reliability", stop_on_failure = TRUE)'
Rscript -e 'devtools::test(filter = "dr-gpsm-errors", stop_on_failure = TRUE)'
```

## Full Test Suite

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
```

Optional learner tests use `testthat::skip_if_not_installed()`. A full
`dr_gpsm()` integration run is performed for each learner available in the
active R library; unavailable learners are reported as skipped.

## Coverage

```sh
Rscript -e 'coverage <- covr::package_coverage(type = "tests"); print(coverage)'
Rscript -e 'covr::report(covr::package_coverage(type = "tests"))'
```

Install the development-only coverage tool with `install.packages("covr")` if
needed.

## Source Package Check

Build first, then check the resulting source archive rather than the working
directory:

```sh
R CMD build .
R CMD check --as-cran BC.GPSM_0.0.1.tar.gz
```

The GitHub Actions workflow runs the source-package check with release R on
macOS, Windows, and Linux.
