# Discrete integration

Computes the discrete integral of a function (or vector of function
evaluations) of a given order, with respect to given design points, and
evaluated at a given query point.

## Usage

``` r
discrete_integ(f, k, xd, x)
```

## Arguments

- f:

  Function, or vector of function evaluations at `c(xd, x)`, the design
  points `xd` adjoined with the query point(s) `x`.

- k:

  Order for the discrete integral calculation. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- x:

  Query point(s).

## Value

Discrete integral of `f` of order `k`, with respect to design points
`xd`, evaluated at the query point(s) `x`.

## Details

The discrete integral operator of order \\k\\, with respect to design
points \\x_1 \< \ldots \< x_n\\, is denoted \\S^k_n\\. It is the inverse
operator to the discrete derivative operator \\\Delta^k_n\\, so that:
\$\$ S^k_n \Delta^k_n = \Delta^k_n S^k_n = \mathrm{Id}, \$\$ where
\\\mathrm{Id}\\ denotes the identity operator. It can also be
represented in a more explicit form, as follows. Acting on a function
\\f\\ of order \\k\\, and evaluated at a query point \\x\\, it yields:
\$\$ (S^k_n f)(x) = \begin{cases} \displaystyle \sum\_{j=1}^k
h^{k-1}\_j(x) \cdot f(x_j) + \sum\_{j=k+1}^i h^{k-1}\_j(x) \cdot
\frac{x_j-x\_{j-k}}{k} \cdot f(x_j) + h^{k-1}\_{i+1}(x) \cdot
\frac{x-x\_{i-k+1}}{k} \cdot f(x) & \\ & \hspace{-75pt} \text{if \$x \in
(x_i,x\_{i+1}\]\$, \$i \geq k\$} \\ \displaystyle \sum\_{j=1}^i
h^{k-1}\_j(x) \cdot f(x_j) \\+\\ h^{k-1}\_{i+1}(x) \cdot f(x) &
\hspace{-75pt} \text{if \$x \in (x_i,x\_{i+1}\]\$, \$i \< k\$} \\ f(x) &
\hspace{-75pt} \text{if \$x \leq x_1\$}, \end{cases} \$\$ where
\\h^{k-1}\_1, \ldots, h^{k-1}\_n\\ denote the falling factorial basis
functions of degree \\k-1\\, with knots in \\x\_{k:(n-1)}\\. The help
file for
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md) gives a
definition of the falling factorial basis. It can be seen (due to the
one-sided support of the falling factorial basis functions) that
discrete integration at \\x = x_i\\, \\i = 1,\ldots,n\\ is equivalent to
multiplication by a weighted version of the falling factorial basis
matrix. For more details, including an alternative recursive
representation for discrete integration (that elucidates its
relationship to discrete differentiation), see Section 3.2 of Tibshirani
(2020).

**Note:** for calculating discrete integrals at the design points
themselves, which could be achieved by taking `x = xd` in the current
function, one should instead use
[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md)
with `di_weighting = TRUE`, as this will be **much** more efficient
(quadratic-time versus linear-time).

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 3.2.

## See also

[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md)
for multiplication by the falling factorial basis matrix, giving a
weighted analog of discrete integration at the design points.

## Examples

``` r
xd = 1:10 / 10
discrete_integ(function(x) 1, 1, xd, xd)
#>  [1] 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9
```
