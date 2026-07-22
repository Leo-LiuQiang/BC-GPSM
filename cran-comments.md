## Submission

This is the first CRAN submission of BC.GPSM.

## Test environments

* Local: macOS 26.5.2, R 4.6.1, arm64
* Win-builder: Windows Server 2022, R-devel (2026-07-20 r90283)
* GitHub Actions: macOS, Windows, and Ubuntu, R release

## R CMD check results

Local `R CMD check --as-cran` results:

* 0 errors
* 0 warnings
* 2 notes

The notes were:

1. This is a new submission.
2. HTML validation was skipped because the external HTML Tidy installed on the
   local system was not recent enough. The HTML manual was validated
   successfully by Win-builder using R-devel.

All 169 test assertions pass locally with no failures, warnings, or skips.
Optional ranger and XGBoost integrations are exercised when those suggested
packages are installed.

The package spelling audit reports no unrecognized words after accounting for
documented statistical abbreviations, model names, and code identifiers. All
URLs in the package metadata and documentation were verified successfully.

## Downstream dependencies

There are currently no downstream CRAN dependencies because this is a new
submission.
