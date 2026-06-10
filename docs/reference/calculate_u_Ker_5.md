# Construct Kernel Matrices up to 5th Order

Extends the recursive kernel construction to fifth order based on four
input kernel matrices. This function builds upon
[`calculate_u_Ker_4()`](https://cxy0714.github.io/HOIF/reference/calculate_u_Ker_4.md).

## Usage

``` r
calculate_u_Ker_5(A1, A2, A3, A4)
```

## Arguments

- A1, A2, A3, A4:

  Numeric square matrices of the same dimension.

## Value

A list containing kernel matrices `Ker2`, `Ker3`, `Ker4`, and `Ker5`.
