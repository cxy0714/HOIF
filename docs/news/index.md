# Changelog

## HOIF 0.2.0

First CRAN submission.

### Adaptation to ustats (\>= 0.1.4)

- The Python dependencies (u-stats, numpy, torch) of the ‘ustats’
  backend are now provisioned automatically on first use (via reticulate
  \>= 1.41 and
  [`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html)
  declared by ‘ustats’), so most users need no manual setup. The
  installation documentation (README, vignette) was rewritten
  accordingly, covering the three setup paths: automatic,
  [`ustats::setup_ustats()`](https://rdrr.io/pkg/ustats/man/setup_ustats.html)
  (CPU-only PyTorch by default, `gpu = TRUE` for CUDA), and
  bring-your-own-environment.
- [`hoif_ate()`](https://cxy0714.github.io/HOIF/reference/hoif_ate.md)
  and
  [`compute_hoif_estimators()`](https://cxy0714.github.io/HOIF/reference/compute_hoif_estimators.md)
  gain a `dtype` argument (“float32”/“float64”, default `NULL` =
  automatic) that is passed through to
  [`ustats::ustat()`](https://rdrr.io/pkg/ustats/man/ustat.html) to
  control the numeric precision of the Python backend.
- The `Remotes:` field and the package-level `Config/reticulate` block
  were removed: ‘ustats’ is now obtained from CRAN and declares its own
  Python requirements.

### Bug fixes

- [`transform_covariates()`](https://cxy0714.github.io/HOIF/reference/transform_covariates.md)
  no longer fails with “object ‘k’ not found” for `method = "splines"`
  and `method = "fourier"` (the basis dimension argument was renamed to
  `basis_dim` without updating the function body). It now also validates
  `basis_dim` and gives an informative error when it is missing or too
  small.
- The Moore-Penrose fallback in
  [`compute_gram_inverse()`](https://cxy0714.github.io/HOIF/reference/compute_gram_inverse.md)
  now checks that ‘MASS’ is installed (it is a suggested package) and
  gives an informative error if not.
- [`print.hoif_ate()`](https://cxy0714.github.io/HOIF/reference/print.hoif_ate.md)
  returns its argument invisibly, as is customary for print methods.

### Documentation

- All exported functions now have executable examples; examples and
  tests that need Python are guarded so that the package can be checked
  without a Python runtime (the pure R backend is exercised instead).
- Added a `tests/testthat.R` runner (previously the test suite was not
  run by `R CMD check`) plus new pure R regression tests.
- The vignette was moved to the standard `vignettes/` directory, its
  code chunks are no longer evaluated at build time (the Python backend
  is not available on build machines), and the references were
  corrected.
- Startup message points users to
  [`ustats::check_ustats_setup()`](https://rdrr.io/pkg/ustats/man/check_ustats_setup.html)
  and the `pure_R_code = TRUE` fallback.
