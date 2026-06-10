# Generate All Ordered m-Tuples of Distinct Indices (Brute Force)

Internal helper function used for brute-force validation of higher-order
influence function (HOIF) calculations. It generates all ordered
\\m\\-permutations of a given index set.

## Usage

``` r
generate_permutations(indices, m)
```

## Arguments

- indices:

  Integer vector of indices.

- m:

  Integer order of the permutation.

## Value

A matrix where each row is an ordered \\m\\-tuple of distinct indices.

## Details

\*\*Warning:\*\* This function has factorial computational complexity
and should only be used for very small sample sizes as part of debugging
or verification routines.
