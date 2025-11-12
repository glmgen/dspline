# In-place computations

Each "dot" function accepts arguments as in its "non-dot" counterpart,
but peforms computations in place, overwriting the first input argument
(which must be a vector of the appropriate length) with the desired
output.

## Usage

``` r
.divided_diff(f, z)

.b_mat_mult(v, k, xd, tf_weighting, transpose, inverse)

.h_mat_mult(v, k, xd, di_weighting, transpose, inverse)
```

## Arguments

- f:

  Function, or vector of function evaluations at the centers.

- z:

  Centers for the divided difference calculation.

- v:

  Vector to be multiplied by B, the extended discrete derivative matrix.

- k:

  Order for the extended discrete derivative matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- tf_weighting:

  Should "trend filtering weighting" be used? This is a weighting of the
  discrete derivatives that is implicit in trend filtering; see details
  for more information. The default is `FALSE`.

- transpose:

  Multiply by the transpose of B? The default is `FALSE`.

- inverse:

  Multiply by the inverse of B? The default is `FALSE`.

- di_weighting:

  Should "discrete integration weighting" be used? Multiplication by
  such a weighted H gives discrete integrals at the design points; see
  details for more information. The default is `FALSE`.

## Value

None. These functions *overwrite* their input.

## Details

These functions should not be used unless you are intentionally doing so
for memory considerations and are nonetheless being very careful.

An **important warning:** each "dot" function only works as expected if
its first argument is passed in as a vector of numeric type. If the
first argument is passed in as an integer vector, then since the output
must (in general) be a numeric vector, it cannot be computed in place
(Rcpp performs an implicit cast and copy when it converts this to
NumericVector type for use in C++).

Also, each "dot" function does not perform any error checking on its
input arguments. Use with care. More details on the computations
performed by individual functions are provided below.

## `.divided_diff()`

Overwrites `f` with all lower-order divided differences: each element
`f[i]` becomes the divided difference with respect to centers `z[1:i]`.
See also
[`divided_diff()`](https://glmgen.github.io/dspline/reference/divided_diff.md).

## `.b_mat_mult()`

Overwrites `v` with `B %*% v`, where `B` is the extended discrete
derivative matrix as returned by
[`b_mat()`](https://glmgen.github.io/dspline/reference/b_mat.md). See
also
[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md).

## `.h_mat_mult()`

Overwrites `v` with `H %*% v`, where `H` is the falling factorial basis
matrix as returned by
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md). See
also
[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md).

## Examples

``` r
v = as.numeric(1:10) # Note: must be of numeric type
b_mat_mult(v, 1, 1:10) 
#>  [1] 1 1 1 1 1 1 1 1 1 1
v
#>  [1]  1  2  3  4  5  6  7  8  9 10
.b_mat_mult(v, 1, 1:10, FALSE, FALSE, FALSE)
v
#>  [1] 1 1 1 1 1 1 1 1 1 1
```
