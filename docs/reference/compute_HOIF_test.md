# Brute-Force Higher-Order Influence Function Estimator (Test Version)

Computes the order-`m` higher-order influence function (HOIF) term using
an explicit enumeration of all ordered index tuples.

## Usage

``` r
compute_HOIF_test(
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

  Binary treatment indicator vector.

- Y:

  Outcome vector.

- mu1:

  Estimated outcome regression under treatment.

- mu0:

  Estimated outcome regression under control.

- pi:

  Estimated propensity scores.

- m:

  Order of the HOIF term.

- sample_splitting:

  Whether to use K-fold sample splitting (0 = no).

- n_folds:

  Number of folds for sample splitting.

- seed:

  Random seed for reproducibility.

## Value

A list with elements:

- IIFF_1_m:

  Order-`m` influence function term for treated units

- IIFF_0_m:

  Order-`m` influence function term for control units

## Details

\*\*This implementation is intentionally naive and is provided solely
for validation and debugging purposes.\*\* It should only be used on
very small datasets to verify the correctness of optimized
implementations.

The computational complexity grows factorially with both sample size and
order `m`.

This function directly implements the combinatorial definition of the
HOIF using full enumeration of index permutations. Matrix inverses are
computed explicitly and no numerical stabilization is included.

This function is part of the package's \*\*internal test
infrastructure\*\* and should not be exported or used in applied
analysis.
