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
#>  $ Omega1: num [1:3, 1:3] 2.257 -0.477 0.362 -0.477 2.333 ...
#>  $ Omega0: num [1:3, 1:3] 1.92 -0.142 -0.197 -0.142 1.573 ...
```
