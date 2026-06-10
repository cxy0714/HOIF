# Build E_j tensor structure

Constructs a nested list representing the tensor structure \[1, \[1,2\],
..., \[j-1,j\], j\] used in Higher-Order Influence Function (HOIF)
calculations.

## Usage

``` r
build_Ej(j)
```

## Arguments

- j:

  Integer greater than or equal to 2 specifying the dimension

## Value

A nested list with j+1 elements representing the tensor structure: -
First element: scalar 1 - Middle elements: vectors \[1,2\], \[2,3\],
..., \[j-1,j\] - Last element: scalar j

## Examples

``` r
build_Ej(3)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 1 2
#> 
#> [[3]]
#> [1] 2 3
#> 
#> [[4]]
#> [1] 3
#> 
```
