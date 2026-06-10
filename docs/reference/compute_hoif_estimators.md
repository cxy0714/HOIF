# Compute HOIF Estimators for ATE

Compute HOIF Estimators for ATE

## Usage

``` r
compute_hoif_estimators(
  residuals,
  B_matrices,
  m = 7,
  backend = "torch",
  pure_R_code = FALSE,
  dtype = NULL
)
```

## Arguments

- residuals:

  A list containing the computed residuals: \`R1\`, \`r1\`, \`R0\`, and
  \`r0\`.

- B_matrices:

  A list containing the projection-like basis matrices: \`B1\` and
  \`B0\`.

- m:

  Integer. The maximum order of the HOIF estimator.

- backend:

  Character. The computation backend used by
  [`ustat`](https://rdrr.io/pkg/ustats/man/ustat.html); either "torch"
  (default) or "numpy". If PyTorch is not available, 'ustats' falls back
  to "numpy" with a warning.

- pure_R_code:

  Logical. Whether to use a native R implementation. This serves as a
  fallback when the Python environment used by 'ustats' (via
  'reticulate') is not available. Note: The pure R implementation only
  supports up to the 6th order (m = 6).

- dtype:

  Optional character string passed to
  [`ustat`](https://rdrr.io/pkg/ustats/man/ustat.html) controlling the
  numeric precision of the Python backend: "float32" or "float64". The
  default `NULL` selects the precision automatically (float32 on a CUDA
  GPU, float64 otherwise). Ignored when `pure_R_code = TRUE`.

## Value

A list of HOIF estimators (ATE, HOIF, IIFF) for orders l = 2, ..., m.

## Examples

``` r
# Pure R example (no Python required), up to order 6
n <- 100
Z <- cbind(1, rnorm(n))
A <- rbinom(n, 1, 0.5)
Y <- A + Z[, 2] + rnorm(n)
residuals <- compute_residuals(A, Y, mu1 = 1 + Z[, 2], mu0 = Z[, 2],
                               pi = rep(0.5, n))
Omega <- compute_gram_inverse(Z, A)
B <- compute_basis_matrix(Z, A, Omega$Omega1, Omega$Omega0)
est <- compute_hoif_estimators(residuals, B, m = 6, pure_R_code = TRUE)
est$ATE
#> [1] -0.04885023 -0.05408448 -0.05178525 -0.05026342 -0.04996561

if (FALSE) { # \dontrun{
# Python backend (requires the 'ustats' Python dependencies), any order
est <- compute_hoif_estimators(residuals, B, m = 7, backend = "torch")
} # }
```
