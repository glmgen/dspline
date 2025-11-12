# Discrete differentiation

Computes the discrete derivative of a function (or vector of function
evaluations) of a given order, with respect to given design points, and
evaluated at a given query point.

## Usage

``` r
discrete_deriv(f, k, xd, x)
```

## Arguments

- f:

  Function, or vector of function evaluations at `c(xd, x)`, the design
  points `xd` adjoined with the query point(s) `x`.

- k:

  Order for the discrete derivative calculation. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- x:

  Query point(s).

## Value

Discrete derivative of `f` of order `k`, with respect to design points
`xd`, evaluated at the query point(s) `x`.

## Details

The discrete derivative operator of order \\k\\, with respect design
points \\x_1 \< \ldots \< x_n\\, is denoted \\\Delta^k_n\\. Acting on a
function \\f\\, and evaluated at a query point \\x\\, it yields: \$\$
(\Delta^k_n f) (x) = \begin{cases} k! \cdot f\[x\_{i-k+1},\ldots,x_i,x\]
& \text{if \$x \in (x_i,x\_{i+1}\]\$, \$i \geq k\$} \\ i! \cdot
f\[x_1,\ldots,x_i,x\] & \text{if \$x \in (x_i,x\_{i+1}\]\$, \$i \< k\$}
\\ f(x) & \text{if \$x \leq x_1\$}, \end{cases} \$\$ where we take
\\x\_{n+1} = \infty\\ for convenience. In other words, for "most" points
\\x \> x_k\\, we define \\(\Delta^k_n f)(x)\\ in terms of a (scaled)
divided difference of \\f\\ of order \\k\\, where the centers are the
\\k\\ points immediately to the left of \\x\\, and \\x\\ itself.
Meanwhile, for "boundary" points \\x \leq x_k\\, we define \\(\Delta^k_n
f)(x)\\ to be a (scaled) divided difference of \\f\\ of the highest
possible order, where the centers are the points to the left of \\x\\,
and \\x\\ itself. For more discussion, including alternative
representations for the discrete differentiation, see Section 3.1 of
Tibshirani (2020).

**Note:** for calculating discrete derivatives at the design points
themselves, which could be achieved by taking `x = xd` in the current
function, one should instead use
[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md)
or
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md),
as these will be more efficient (both will be linear-time, but the
latter functions will be faster).

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 3.1.

## See also

[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md),
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md)
for multiplication by the extended and non-extended discrete derivative
matrices, giving discrete derivatives at design points.

## Examples

``` r
xd = 1:10 / 10
discrete_deriv(function(x) x^2, 1, xd, xd)
#>  [1] 0.01 0.30 0.50 0.70 0.90 1.10 1.30 1.50 1.70 1.90
```
