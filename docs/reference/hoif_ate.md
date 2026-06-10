# Main function: HOIF estimators for ATE with optional sample splitting

Computes the higher-order influence function terms of orders 2 to m,
which estimate the estimable bias of the standard first-order doubly
robust (AIPW) estimator of the ATE and are used to debias it.

## Usage

``` r
hoif_ate(
  X,
  A,
  Y,
  mu1,
  mu0,
  pi,
  transform_method = "none",
  basis_dim = NULL,
  inverse_method = "direct",
  m = 7,
  sample_split = FALSE,
  n_folds = 2,
  backend = "torch",
  seed = 42,
  pure_R_code = FALSE,
  dtype = NULL,
  ...
)
```

## Arguments

- X:

  Covariate matrix (n x p)

- A:

  Treatment vector (n x 1)

- Y:

  Outcome vector (n x 1)

- mu1:

  Predicted outcomes under treatment (predictions supplied by the user,
  ideally estimated on a separate, independent sample)

- mu0:

  Predicted outcomes under control (see \`mu1\`)

- pi:

  Predicted propensity scores (see \`mu1\`)

- transform_method:

  Character: method to transform covariates before constructing basis
  functions. - "splines": use basis splines expansion - "fourier": use
  Fourier basis expansion - "none": no transformation (use raw
  covariates)

- basis_dim:

  Integer: number of basis functions to generate when using "splines" or
  "fourier" transformations. Higher values provide more flexible
  approximations but may increase variance.

- inverse_method:

  Character: regularization method for Gram matrix inversion. -
  "direct": Cholesky-based inversion (falls back to the Moore-Penrose
  inverse from MASS when the Gram matrix is not positive definite) -
  "nlshrink": nonlinear shrinkage estimator (Ledoit-Wolf type) -
  "corpcor": shrinkage via the corpcor package (for high-dimensional
  settings)

- m:

  Maximum order for HOIF (up to 6 when `pure_R_code = TRUE`)

- sample_split:

  Logical: whether to cross-fit the inverse weighted Gram matrix against
  the U-statistics. If \`TRUE\` (the eHOIF case), the sample is split
  into \`n_folds\` folds; for each fold, the Gram matrix is estimated on
  the remaining folds, the U-statistics are computed on that fold, and
  the results are averaged across folds. If \`FALSE\` (the sHOIF case),
  both are computed on the same sample, without distinction. Note this
  is not a cross-fitting of the nuisance functions, whose predictions
  are supplied via \`mu1\`, \`mu0\`, \`pi\`.

- n_folds:

  Number of folds for sample splitting (if used)

- backend:

  Character: computation backend used by
  [`ustat`](https://rdrr.io/pkg/ustats/man/ustat.html); "torch"
  (default) or "numpy". Ignored when `pure_R_code = TRUE`.

- seed:

  Random seed for reproducibility (for sample splitting)

- pure_R_code:

  Logical: if \`TRUE\`, the higher-order U-statistics are computed with
  a pure R implementation (no Python required), which supports orders up
  to m = 6. If \`FALSE\` (default), they are computed by the 'ustats'
  package, whose Python dependencies are provisioned automatically on
  first use (see the package README).

- dtype:

  Optional character string ("float32" or "float64") controlling the
  numeric precision of the Python backend; \`NULL\` (default) selects
  the precision automatically. Passed to
  [`ustat`](https://rdrr.io/pkg/ustats/man/ustat.html); ignored when
  `pure_R_code = TRUE`.

- ...:

  Additional arguments passed to transform_covariates

## Value

An object of class `"hoif_ate"`: a list with components

- ATE:

  ATE estimates for orders 2 to m

- HOIF1, HOIF0:

  HOIF estimates for the treated/control arm

- IIFF1, IIFF0:

  Incremental influence function terms for the treated/control arm

- orders:

  The orders 2 to m

- convergence_data:

  Data frame with the ATE estimate per order

## Details

Conceptually, HOIF estimation involves three estimation tasks, and
ideally each uses its own, independent part of the data: (1) estimating
the nuisance functions mu(1, X), mu(0, X) and pi(X); (2) estimating the
inverse weighted Gram matrix; (3) computing the higher-order
U-statistics. This package does not implement task (1): \`hoif_ate()\`
only takes the nuisance \*predictions\* \`mu1\`, \`mu0\` and \`pi\` as
inputs, so the overall three-way cross-fitting is left to the user. The
\`sample_split\` argument controls only the split between tasks (2) and
(3); see its description below.

## See also

[`compute_hoif_estimators`](https://cxy0714.github.io/HOIF/reference/compute_hoif_estimators.md)
for the lower-level estimation routine.

## Examples

``` r
# A small, self-contained example using the pure R backend
set.seed(1)
n <- 100
X <- matrix(rnorm(n), ncol = 1)
A <- rbinom(n, 1, 0.5)
Y <- as.numeric(A + X %*% 0.5 + rnorm(n, 0, 0.1))
mu1 <- as.numeric(1 + X %*% 0.5)
mu0 <- as.numeric(X %*% 0.5)
pi <- rep(0.5, n)

fit <- hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
                transform_method = "none", m = 6,
                pure_R_code = TRUE)
print(fit)
#> HOIF Estimators for Average Treatment Effect
#> =============================================
#> 
#> Higher-order correction terms by order:
#>   Order    ATE  HOIF1 HOIF0
#> 1     2 -5e-04 -3e-04 1e-04
#> 2     3 -4e-04 -3e-04 2e-04
#> 3     4 -4e-04 -2e-04 2e-04
#> 4     5 -4e-04 -2e-04 2e-04
#> 5     6 -4e-04 -2e-04 2e-04
#> 
#> Estimated AIPW bias correction for the ATE (highest order): -4e-04 
#> (add this value to the first-order AIPW/DR estimate of the ATE to debias it)

if (FALSE) { # \dontrun{
# Python backend (provisioned automatically on first use), order m = 7,
# with 2-fold sample splitting
fit <- hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
                transform_method = "none", m = 7,
                sample_split = TRUE, n_folds = 2, seed = 123,
                backend = "torch")
print(fit)
plot(fit)
} # }
```
