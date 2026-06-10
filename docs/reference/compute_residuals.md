# Compute residuals for both treatment groups

Compute residuals for both treatment groups

## Usage

``` r
compute_residuals(A, Y, mu1, mu0, pi)
```

## Arguments

- A:

  Treatment vector (0 or 1)

- Y:

  Outcome vector

- mu1:

  Predicted outcomes under treatment (mu(1, X))

- mu0:

  Predicted outcomes under control (mu(0, X))

- pi:

  Propensity scores

## Value

List with R1, r1, R0, r0

## Examples

``` r
n <- 100
A <- rbinom(n, 1, 0.5)
Y <- rnorm(n)
res <- compute_residuals(A, Y, mu1 = rep(0.5, n), mu0 = rep(-0.5, n),
                         pi = rep(0.5, n))
str(res)
#> List of 4
#>  $ R1: num [1:100] 1.035 -0.916 -1.021 0.351 -0.166 ...
#>  $ r1: num [1:100] -1 -1 1 -1 -1 1 1 1 -1 1 ...
#>  $ R0: num [1:100] 2.0349 0.0838 -0.0205 1.3506 0.8345 ...
#>  $ r0: num [1:100] 1 1 -1 1 1 -1 -1 -1 1 -1 ...
```
