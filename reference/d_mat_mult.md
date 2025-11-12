# Multiply by D matrix

Multiplies a given vector by D, the discrete derivative matrix of a
given order, with respect to given design points.

## Usage

``` r
d_mat_mult(v, k, xd, tf_weighting = FALSE, transpose = FALSE)
```

## Arguments

- v:

  Vector to be multiplied by D, the discrete derivative matrix.

- k:

  Order for the discrete derivative matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- tf_weighting:

  Should "trend filtering weighting" be used? This is a weighting of the
  discrete derivatives that is implicit in trend filtering; see details
  for more information. The default is `FALSE`.

- transpose:

  Multiply by the transpose of D? The default is `FALSE`.

## Value

Product of the discrete derivative matrix D and the input vector `v`.

## Details

The discrete derivative matrix of order \\k\\, with respect to design
points \\x_1 \< \ldots \< x_n\\, is denoted \\D^k_n\\. It has dimension
\\(n-k) \times n\\. Acting on a vector \\v\\ of function evaluations at
the design points, denoted \\v = f(x\_{1:n})\\, it gives the discrete
derivatives of \\f\\ at the points \\x\_{(k+1):n}\\: \$\$ D^k_n v =
(\Delta^k_n f) (x\_{(k+1):n}). \$\$ The matrix \\D^k_n\\ can be
constructed recursively as the product of a diagonally-weighted first
difference matrix and \\D^{k-1}\_n\\; see the help file for
[`d_mat()`](https://glmgen.github.io/dspline/reference/d_mat.md), or
Section 6.1 of Tibshirani (2020). Therefore, multiplication by \\D^k_n\\
or by its transpose can be performed in \\O(nk)\\ operations based on
iterated weighted differences. See Appendix D of Tibshirani (2020) for
details.

The option `tf_weighting = TRUE` performs multiplication by \\W^k_n
D^k_n\\ where \\W^k_n\\ is a \\(n-k) \times (n-k)\\ diagonal matrix with
entries \\(x\_{i+k} - x_i) / k\\, \\i = 1,\ldots,n-k\\. This weighting
is implicit in trend filtering, as the penalty in the \\k\\th order
trend filtering optimization problem (with optimization parameter
\\\theta\\) is \\\\W^{k+1}\_n D^{k+1}\_n \theta\\\_1\\. Moreover, this
is precisely the \\k\\th order total variation of the \\k\\th degree
discrete spline interpolant \\f\\ to \\\theta\\, with knots in
\\x\_{(k+1):(n-1)}\\; that is, such an interpolant satisfies: \$\$
\mathrm{TV}(D^k f) = \\W^{k+1}\_n D^{k+1}\_n \theta\\\_1, \$\$ where
\\D^k f\\ is the \\k\\th derivative of \\f\\. See Section 9.1. of
Tibshirani (2020) for more details.

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.1.

## See also

[`discrete_deriv()`](https://glmgen.github.io/dspline/reference/discrete_deriv.md)
for discrete differentiation at arbitrary query points,
[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md)
for multiplying by the extended discrete derivative matrix, and
[`d_mat()`](https://glmgen.github.io/dspline/reference/d_mat.md) for
constructing the discrete derivative matrix.

## Examples

``` r
v = sort(runif(10))
as.vector(d_mat(2, 1:10) %*% v)
#> [1] -0.050574895  0.004873456  0.042702133  0.006656813 -0.091228792
#> [6] -0.050737241  0.227824376 -0.229034467
d_mat_mult(v, 2, 1:10) 
#> [1] -0.050574895  0.004873456  0.042702133  0.006656813 -0.091228792
#> [6] -0.050737241  0.227824376 -0.229034467
```
