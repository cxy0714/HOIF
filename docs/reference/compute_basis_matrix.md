# Compute basis projection matrices

Compute basis projection matrices

## Usage

``` r
compute_basis_matrix(Z, A, Omega1, Omega0)
```

## Arguments

- Z:

  Basis matrix (n x k)

- A:

  Treatment vector (n x 1)

- Omega1:

  Inverse Gram matrix for treatment group

- Omega0:

  Inverse Gram matrix for control group

## Value

List with B1 and B0 (projection matrices)

## Examples

``` r
n <- 100
Z <- cbind(1, matrix(rnorm(n * 2), n, 2))
A <- rbinom(n, 1, 0.5)
Omega <- compute_gram_inverse(Z, A)
B <- compute_basis_matrix(Z, A, Omega$Omega1, Omega$Omega0)
dim(B$B1)
#> [1] 100 100
```
