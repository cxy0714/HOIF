## New submission (0.2.0)

This is the first submission of the HOIF package to CRAN.

HOIF implements Higher-Order Influence Function estimators of the Average
Treatment Effect. The higher-order U-statistics at its core are computed by
the 'ustats' package (an R interface to the Python package 'u-stats'); a
pure R implementation (orders up to 6) is provided as a fallback that does
not require Python.

Note: HOIF imports 'ustats' (>= 0.1.4), which is currently under review at
CRAN. This submission depends on 'ustats' being accepted first.

## Test environments

* Local Linux (Ubuntu 24.04, R 4.3.3)
* (Please update with rhub / win-builder results before submitting)

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: New submission.

## Comments on examples, tests and vignettes

* Examples that require the Python runtime are wrapped in \dontrun{},
  because the Python dependencies ('u-stats', 'numpy', 'torch') are not
  available on CRAN check machines and would trigger a large one-time
  download. Every exported function additionally has an executable example
  that runs without Python (using the pure R backend).
* For the same reason, Python-dependent tests call skip_on_cran() in
  addition to runtime availability checks, so the test suite never
  initializes Python on CRAN machines. The pure R backend (including a
  brute-force correctness check of the estimator) is tested
  unconditionally.
* The vignette sets eval = FALSE globally because most of its code requires
  the Python runtime.

## Downstream dependencies

None. This is a new package.
