# HOIF Estimators for Average Treatment Effect

This file implements the Higher Order Influence Function (HOIF)
estimators for Average Treatment Effect (ATE) estimation with nuisance
functions.

## Usage

``` r
transform_covariates(X, method = "splines", basis_dim, degree = 3, period = 1)
```

## Arguments

- X:

  Matrix of covariates (n x p)

- method:

  Character: "splines", "fourier", or "none"

- basis_dim:

  Integer: dimension of basis expansion per covariate (ignored if method
  = "none"; must be at least degree + 1 for "splines" and at least 2 for
  "fourier")

- degree:

  Integer: degree for B-splines (default 3; ignored if method !=
  "splines")

- period:

  Numeric: period for Fourier basis (default 1; ignored if method !=
  "fourier")

## Value

Matrix of transformed covariates. For method = "none" the input X is
returned unchanged (n x p); for "splines" and "fourier" the basis
expansions of all covariates are column-bound together with an intercept
column.

## Author

Xingyu Chen Transform covariates to basis functions

## Examples

``` r
X <- matrix(rnorm(40), nrow = 20, ncol = 2)
Z_splines <- transform_covariates(X, method = "splines", basis_dim = 5)
Z_fourier <- transform_covariates(X, method = "fourier", basis_dim = 4)
Z_raw <- transform_covariates(X, method = "none")
dim(Z_splines)
#> [1] 20  9
dim(Z_fourier)
#> [1] 20  9
```
