# Construct N matrix

Constructs the discrete B-spline basis matrix of a given order, with
respect to given design points and given knot points.

## Usage

``` r
n_mat(k, xd, normalized = TRUE, knot_idx = NULL)
```

## Arguments

- k:

  Order for the discrete B-spline basis matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- normalized:

  Should the discrete B-spline basis vectors be normalized to attain a
  maximum value of 1 over the design points? The default is `TRUE`.

- knot_idx:

  Vector of indices, a subset of `(k+1):(n-1)` where `n = length(xd)`,
  that indicates which design points should be used as knot points for
  the discrete B-splines. Must be sorted in increasing order. The
  default is `NULL`, which is taken to mean `(k+1):(n-1)` as in the
  "canonical" discrete spline space (though in this case the returned N
  matrix will be trivial: it will be the identity matrix). See details.

## Value

Sparse matrix of dimension `length(xd)` by `length(knot_idx) + k + 1`.

## Details

The discrete B-spline basis matrix of order \\k\\, with respect to
design points \\x_1 \< \ldots \< x_n\\, and knot set \\T \subseteq
x\_{(k+1):(n-1)}\\ is denoted \\N^k_T\\. It has dimension \\(\|T\|+k+1)
\times n\\, and its entries are given by: \$\$ (N^k_T)\_{ij} =
\eta^k_j(x_i), \$\$ where \\\eta^k_1, \ldots, \eta^k_m\\ are discrete
B-spline (DB-spline) basis functions and \\m = \|T\|+k+1\\. As is
suggested by their name, the DB-spline functions are linearly
independent and span the space of discrete splines with knots at \\T\\.
Each DB-spline \\\eta^k_j\\ has a key local support property: it is
supported on an interval containing at most \\k+2\\ adjacent knots.

The functions \\\eta^k_1, \ldots, \eta^k_m\\ are, in general, not
available in closed-form, and are defined by setting up and solving a
sequence of locally-defined linear systems. For any knot set \\T\\,
computation of the evaluations of all DB-splines at the design points
can be done in \\O(nk^2)\\ operations; see Sections 7, 8.2, and 8.3 of
Tibshirani (2020) for details. The current function uses a sparse QR
decomposition from the `Eigen::SparseQR` module in C++ in order to solve
the local linear systems.

When \\T = x\_{(k+1):(n-1)}\\, the knot set corresponding to the
"canonical" discrete spline space (spanned by the falling factorial
basis functions \\h^k_1, \ldots, h^k_n\\ whose evaluations make up
\\H^k_n\\; see the help file for
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md)), the
DB-spline basis matrix, which we denote by \\N^k_n\\, is trivial: it
equals the \\n \times n\\ identity matrix, \\N^k_n = I_n\\. Therefore
DB-splines are really only useful for knot sets \\T\\ that are proper
subsets of \\x\_{(k+1):(n-1)}\\. Specification of the knot set \\T\\ is
done via the argument `knot_idx`.

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Sections 7, 8.2, and 8.3.

## See also

[`h_eval()`](https://glmgen.github.io/dspline/reference/h_eval.md) for
constructing evaluations of the discrete B-spline basis at arbitrary
query points.

## Examples

``` r
n_mat(2, 1:10, knot_idx = c(3, 5, 7))
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
n_mat(2, 1:10, knot_idx = c(4, 6, 8))
#> 10 x 6 sparse Matrix of class "dgCMatrix"
#>                                                                  
#>  [1,] 1.0000000 .         .         .         .         .        
#>  [2,] 0.3333333 0.8888889 .         .         .         .        
#>  [3,] .         1.0000000 0.3333333 .         .         .        
#>  [4,] .         0.3333333 1.0000000 .         .         .        
#>  [5,] .         .         1.0000000 0.3333333 .         .        
#>  [6,] .         .         0.3333333 1.0000000 .         .        
#>  [7,] .         .         .         1.0000000 0.3333333 .        
#>  [8,] .         .         .         0.3333333 1.0000000 .        
#>  [9,] .         .         .         .         0.8888889 0.3333333
#> [10,] .         .         .         .         .         1.0000000
```
