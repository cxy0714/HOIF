# Compute U-statistics HOIF-type from Order 2 to 6 in pure R code.

This function serves as a core computational component in higher-order
influence function (HOIF) estimators in pure R code.

## Usage

``` r
calculate_u_statistics_pure_r_six(Vector_1, Vector_2, A1, A2, A3, A4, A5)
```

## Arguments

- Vector_1:

  Numeric column vector of length \\n\\.

- Vector_2:

  Numeric column vector of length \\n\\.

- A1, A2, A3, A4, A5:

  Numeric \\n \times n\\ kernel matrices of the same dimension.

## Value

A named list containing numeric U-statistic estimates:

- U2:

  Second-order U-statistic

- U3:

  Third-order U-statistic

- U4:

  Fourth-order U-statistic

- U5:

  Fifth-order U-statistic

- U6:

  Sixth-order U-statistic

## Details

Internally, the function constructs kernel matrices for orders 2 through
6 using recursive matrix operations and removes diagonal contributions
to ensure degenerate U-statistics.

All diagonal elements of intermediate kernel matrices are removed to
avoid self-interactions. Matrix multiplications are performed via
\`eigenMapMatMult()\` and element-wise products via \`hadamard()\`. The
exact formula of the output is: \$\$ \mathbb{U}\_{n,m} =
\frac{1}{\binom{n}{m} m!} \sum\_{i_1 \ne \cdots \ne i_m} Vector_1\[i_1\]
\cdot A1\[i_1,i_2\] \cdot A1\[i_2,i_3\] \cdots A1\[i\_{m-1},i\_{m}\]
\cdot Vector_2\[i\_{m}\] \$\$
