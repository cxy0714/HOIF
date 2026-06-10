# Generate All Permutations of a Vector (Recursive, Brute Force)

Internal recursive helper used by
[`generate_permutations()`](https://cxy0714.github.io/HOIF/reference/generate_permutations.md)
to enumerate all permutations of a vector.

## Usage

``` r
permn(x)
```

## Arguments

- x:

  A vector.

## Value

A list of vectors, each being a permutation of `x`.

## Details

\*\*Warning:\*\* Extremely computationally expensive for vectors longer
than ~8 elements. Only intended for internal brute-force validation
code.
