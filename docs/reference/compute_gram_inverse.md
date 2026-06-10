# Compute inverse of weighted Gram matrix

Compute inverse of weighted Gram matrix

## Usage

``` r
compute_gram_inverse(Z, A, method = "direct")
```

## Arguments

- Z:

  Basis matrix (n x k)

- A:

  Treatment vector

- method:

  Character: "direct" (Cholesky-based inversion, falling back to the
  Moore-Penrose inverse from MASS if the Gram matrix is not positive
  definite), "nlshrink" (nonlinear shrinkage), or "corpcor"
  (pseudoinverse via corpcor)

## Value

List with Omega1 and Omega0 (inverse Gram matrices)

## Examples

``` r
n <- 100
Z <- cbind(1, matrix(rnorm(n * 2), n, 2))
A <- rbinom(n, 1, 0.5)
Omega <- compute_gram_inverse(Z, A)
str(Omega)
#> List of 2
#>  $ Omega1: num [1:3, 1:3] 2.165 -0.119 -0.236 -0.119 1.489 ...
#>  $ Omega0: num [1:3, 1:3] 1.96 -0.311 0.267 -0.311 2.278 ...
```
