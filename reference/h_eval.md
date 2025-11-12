# Evaluate H basis

Evaluates the falling factorial basis of a given order, with respect to
given design points, at arbitrary query points.

## Usage

``` r
h_eval(k, xd, x, col_idx = NULL)
```

## Arguments

- k:

  Order for the falling factorial basis. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- x:

  Query points. Must be sorted in increasing order.

- col_idx:

  Vector of indices, a subset of `1:n` where `n = length(xd)`, that
  indicates which columns of the constructed matrix should be returned.
  The default is `NULL`, which is taken to mean `1:n`.

## Value

Sparse matrix of dimension `length(x)` by `length(col_idx)`.

## Details

The falling factorial basis functions of order \\k\\, defined with
respect to design points \\x_1 \< \ldots \< x_n\\, are denoted \\h^k_1,
\ldots, h^k_n\\. For their precise definition and further references,
see the help file for
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md). The
current function produces a matrix of evaluations of the falling
factorial basis at an arbitrary sequence of query points. For each query
point \\x\\, this matrix has a corresponding row with entries: \$\$
h^k_j(x), \\ j = 1, \ldots, n. \$\$

## See also

[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md) for
constructing evaluations of the falling factorial basis at the design
points.

## Examples

``` r
xd = 1:10 / 10
x = 1:9 / 10 + 0.05
h_mat(2, xd)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                    
#>  [1,] 1 .   .    .    .    .    .    .    .    .   
#>  [2,] 1 0.1 .    .    .    .    .    .    .    .   
#>  [3,] 1 0.2 0.01 .    .    .    .    .    .    .   
#>  [4,] 1 0.3 0.03 0.01 .    .    .    .    .    .   
#>  [5,] 1 0.4 0.06 0.03 0.01 .    .    .    .    .   
#>  [6,] 1 0.5 0.10 0.06 0.03 0.01 .    .    .    .   
#>  [7,] 1 0.6 0.15 0.10 0.06 0.03 0.01 .    .    .   
#>  [8,] 1 0.7 0.21 0.15 0.10 0.06 0.03 0.01 .    .   
#>  [9,] 1 0.8 0.28 0.21 0.15 0.10 0.06 0.03 0.01 .   
#> [10,] 1 0.9 0.36 0.28 0.21 0.15 0.10 0.06 0.03 0.01
h_eval(2, xd, x)
#> 9 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                              
#>  [1,] 1 0.05 -0.00125 .       .       .       .       .       .       .      
#>  [2,] 1 0.15  0.00375 .       .       .       .       .       .       .      
#>  [3,] 1 0.25  0.01875 0.00375 .       .       .       .       .       .      
#>  [4,] 1 0.35  0.04375 0.01875 0.00375 .       .       .       .       .      
#>  [5,] 1 0.45  0.07875 0.04375 0.01875 0.00375 .       .       .       .      
#>  [6,] 1 0.55  0.12375 0.07875 0.04375 0.01875 0.00375 .       .       .      
#>  [7,] 1 0.65  0.17875 0.12375 0.07875 0.04375 0.01875 0.00375 .       .      
#>  [8,] 1 0.75  0.24375 0.17875 0.12375 0.07875 0.04375 0.01875 0.00375 .      
#>  [9,] 1 0.85  0.31875 0.24375 0.17875 0.12375 0.07875 0.04375 0.01875 0.00375
```
