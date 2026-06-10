# Brute-Force Sequence of Higher-Order Influence Function Estimates (Test Only)

Computes a sequence of higher-order influence function (HOIF) estimates
from order 2 up to order `m` using a brute-force implementation. This
function is \*\*only intended for double-checking correctness\*\* of the
main fast implementation.

## Usage

``` r
compute_HOIF_sequence_test(
  X,
  A,
  Y,
  mu1,
  mu0,
  pi,
  m,
  sample_splitting = 0,
  n_folds = 2,
  seed = 42
)
```

## Arguments

- X:

  Covariate matrix (n x p).

- A:

  Binary treatment vector of length n.

- Y:

  Outcome vector of length n.

- mu1:

  Estimated outcome regression under treatment.

- mu0:

  Estimated outcome regression under control.

- pi:

  Estimated propensity scores.

- m:

  Maximum order of HOIF to compute.

- sample_splitting:

  Logical or integer; whether to use sample splitting.

- n_folds:

  Number of folds for sample splitting.

- seed:

  Random seed for fold assignment.

## Value

A list containing cumulative HOIF and incremental IIFF terms for both
treatment arms.

## Details

\*\*Warning:\*\* The computation scales combinatorially with both sample
size and order and becomes infeasible beyond very small datasets.

This routine repeatedly calls
[`compute_HOIF_test()`](https://cxy0714.github.io/HOIF/reference/compute_HOIF_test.md)
and accumulates incremental influence function terms. It is not
optimized and should never be used in production or large-sample
simulations.
