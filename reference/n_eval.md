# Evaluate N basis

Evaluates the discrete B-spline basis of a given order, with respect to
given design points, evaluated at arbitrary query points.

## Usage

``` r
n_eval(k, xd, x, normalized = TRUE, knot_idx = NULL, N = NULL)
```

## Arguments

- k:

  Order for the discrete B-spline basis. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- x:

  Query points. Must be sorted in increasing order.

- normalized:

  Should the discrete B-spline basis vectors be normalized to attain a
  maximum value of 1 over the design points? The default is `TRUE`.

- knot_idx:

  Vector of indices, a subset of `(k+1):(n-1)` where `n = length(xd)`,
  that indicates which design points should be used as knot points for
  the discrete B-splines. Must be sorted in increasing order. The
  default is `NULL`, which is taken to mean `(k+1):(n-1)`.

- N:

  Matrix of discrete B-spline evaluations at the design points. The
  default is `NULL`, which means that this is precomputed before
  constructing the matrix of discrete B-spline evaluations at the query
  points. If `N` is non-`NULL`, then the argument `normalized` will be
  ignored (as this would have only been used to construct N at the
  design points).

## Value

Sparse matrix of dimension `length(x)` by `length(knot_idx) + k + 1`.

## Details

The discrete B-spline basis functions of order \\k\\, defined with
respect to design points \\x_1 \< \ldots \< x_n\\, are denoted
\\\eta^k_1, \ldots, \eta^k_n\\. For a discussion of their properties and
further references, see the help file for
[`n_mat()`](https://glmgen.github.io/dspline/reference/n_mat.md). The
current function produces a matrix of evaluations of the discrete
B-spline basis at an arbitrary sequence of query points. For each query
point \\x\\, this matrix has a corresponding row with entries: \$\$
\eta^k_j(x), \\ j = 1, \ldots, n. \$\$

Unlike the falling factorial basis, the discrete B-spline basis is not
generally available in closed-form. Therefore, the current function
(unlike
[`h_eval()`](https://glmgen.github.io/dspline/reference/h_eval.md)) will
first check if it should precompute the evaluations of the discrete
B-spline basis at the design points. If the argument `N` is non-`NULL`,
then it will use this as the matrix of evaluations at the design points;
if `N` is `NULL`, then it will call
[`n_mat()`](https://glmgen.github.io/dspline/reference/n_mat.md) to
produce such a matrix, and will pass to this function the arguments
`normalized` and `knot_idx` accordingly.

After obtaining the matrix of discrete B-spline evaluations at the
design points, the fast interpolation scheme from
[`dspline_interp()`](https://glmgen.github.io/dspline/reference/dspline_interp.md)
is used to produce evaluations at the query points.

## See also

[`n_mat()`](https://glmgen.github.io/dspline/reference/n_mat.md) for
constructing evaluations of the discrete B-spline basis at the design
points.

## Examples

``` r
xd = 1:10 / 10
x = 1:9 / 10 + 0.05
n_mat(2, xd, knot_idx = c(3, 5, 7))
#> 10 x 6 sparse Matrix of class "dgCMatrix"
#>                                                          
#>  [1,] 1 .         .         .         .         .        
#>  [2,] . 1.0000000 .         .         .         .        
#>  [3,] . 0.3333333 0.8888889 .         .         .        
#>  [4,] . .         1.0000000 0.3125000 .         .        
#>  [5,] . .         0.3333333 0.9375000 .         .        
#>  [6,] . .         .         1.0000000 0.2857143 .        
#>  [7,] . .         .         0.5000000 0.8571429 .        
#>  [8,] . .         .         0.1666667 1.0000000 0.1666667
#>  [9,] . .         .         .         0.7142857 0.5000000
#> [10,] . .         .         .         .         1.0000000
n_eval(2, xd, x, knot_idx = c(3, 5, 7))
#> 9 x 6 sparse Matrix of class "dgCMatrix"
#>                                                                     
#>  [1,]  0.375  0.70833333 -0.11111111  .          .         .        
#>  [2,] -0.125  0.87500000  0.33333333  .          .         .        
#>  [3,]  .      0.12500000  1.04166667  0.11718750 .         .        
#>  [4,]  .     -0.04166667  0.76388889  0.58593750 .         .        
#>  [5,]  .      .           0.12500000  1.03906250 0.1071429 .        
#>  [6,]  .      .          -0.04166667  0.82031250 0.5357143 .        
#>  [7,]  .      .           .           0.31250000 0.9821429 0.0625000
#>  [8,]  .      .           .           0.06250000 0.9107143 0.3125000
#>  [9,]  .      .           .          -0.02083333 0.4107143 0.7291667
```
