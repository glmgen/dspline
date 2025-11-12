# Multiply by B matrix

Multiplies a given vector by B, the extended discrete derivative matrix
of a given order, with respect to given design points.

## Usage

``` r
b_mat_mult(v, k, xd, tf_weighting = FALSE, transpose = FALSE, inverse = FALSE)
```

## Arguments

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

## Value

Product of the extended discrete derivative matrix B and the input
vector `v`.

## Details

The extended discrete derivative matrix of order \\k\\, with respect to
design points \\x_1 \< \ldots \< x_n\\, is denoted \\B^k_n\\. It is
square, having dimension \\n \times n\\. Acting on a vector \\v\\ of
function evaluations at the design points, denoted \\v = f(x\_{1:n})\\,
it gives the discrete derivatives of \\f\\ at the points \\x\_{1:n}\\:
\$\$ B^k_n v = (\Delta^k_n f) (x\_{1:n}). \$\$ The matrix \\B^k_n\\ can
be constructed recursively as the product of a diagonally-weighted first
difference matrix and \\B^{k-1}\_n\\; see the help file for
[`b_mat()`](https://glmgen.github.io/dspline/reference/b_mat.md), or
Section 6.2 of Tibshirani (2020). Therefore, multiplication by \\B^k_n\\
or by its transpose can be performed in \\O(nk)\\ operations based on
iterated weighted differences. See Appendix D of Tibshirani (2020) for
details.

The option `tf_weighting = TRUE` performs multiplication by \\Z^k_n
B^k_n\\ where \\Z^k_n\\ is an \\n \times n\\ diagonal matrix whose top
left \\k \times k\\ block equals the identity matrix and bottom right
\\(n-k) \times (n-k)\\ block equals \\W^k_n\\, the latter being a
diagonal weight matrix that is implicit in trend filtering, as explained
in the help file for
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md).

Lastly, the matrix \\B^k_n\\ has a special **inverse relationship** to
the falling factorial basis matrix \\H^{k-1}\_n\\ of degree \\k-1\\ with
knots in \\x\_{k:(n-1)}\\; it satisfies: \$\$ Z^k_n B^k_n H^{k-1}\_n =
I_n, \$\$ where \\Z^k_n\\ is the \\n \times n\\ diagonal matrix as
described above, and \\I_n\\ is the \\n \times n\\ identity matrix.
This, combined with the fact that the falling factorial basis matrix has
an efficient recursive representation in terms of weighted cumulative
sums, means that multiplying by \\(B^k_n)^{-1}\\ or its transpose can be
performed in \\O(nk)\\ operations. See Section 6.3 and Appendix D of
Tibshirani (2020) for details.

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.2.

## See also

[`discrete_deriv()`](https://glmgen.github.io/dspline/reference/discrete_deriv.md)
for discrete differentiation at arbitrary query points,
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md)
for multiplying by the discrete derivative matrix, and
[`b_mat()`](https://glmgen.github.io/dspline/reference/b_mat.md) for
constructing the extended discrete derivative matrix.

## Examples

``` r
v = sort(runif(10))
as.vector(b_mat(2, 1:10) %*% v)
#>  [1]  0.007399441  0.026841891  0.113857403 -0.025872676 -0.084208132
#>  [6]  0.115389280 -0.114623875  0.203720707 -0.195465074  0.062439625
b_mat_mult(v, 2, 1:10) 
#>  [1]  0.007399441  0.026841891  0.113857403 -0.025872676 -0.084208132
#>  [6]  0.115389280 -0.114623875  0.203720707 -0.195465074  0.062439625
```
