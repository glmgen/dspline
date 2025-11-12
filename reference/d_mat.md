# Construct D matrix

Constructs the discrete derivative matrix of a given order, with respect
to given design points.

## Usage

``` r
d_mat(k, xd, tf_weighting = FALSE, row_idx = NULL)
```

## Arguments

- k:

  Order for the discrete derivative matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- tf_weighting:

  Should "trend filtering weighting" be used? This is a weighting of the
  discrete derivatives that is implicit in trend filtering; see details
  for more information. The default is `FALSE`.

- row_idx:

  Vector of indices, a subset of `1:(n-k)` where `n = length(xd)`, that
  indicates which rows of the constructed matrix should be returned. The
  default is `NULL`, which is taken to mean `1:(n-k)`.

## Value

Sparse matrix of dimension `length(row_idx)` by `length(xd)`.

## Details

The discrete derivative matrix of order \\k\\, with respect to design
points \\x_1 \< \ldots \< x_n\\, is denoted \\D^k_n\\. It has dimension
\\(n-k) \times n\\, and is banded with bandwidth \\k+1\\. It can be
constructed recursively, as follows. We first define the \\(n-1) \times
n\\ first difference matrix \\\bar{D}\_n\\: \$\$ \bar{D}\_n =
\left\[\begin{array}{rrrrrr} -1 & 1 & 0 & \ldots & 0 & 0 \\ 0 & -1 & 1 &
\ldots & 0 & 0 \\ \vdots & & & & & \\ 0 & 0 & 0 & \ldots & -1 & 1
\end{array}\right\], \$\$ and for \\k \geq 1\\, define the \\(n-k)
\times (n-k)\\ diagonal weight matrix \\W^k_n\\ to have diagonal entries
\\(x\_{i+k} - x_i) / k\\, \\i = 1,\ldots,n-k\\. The \\k\\th order
discrete derivative matrix \\D^k_n\\ is then given by the recursion:
\$\$ \begin{aligned} D^1_n &= (W^1_n)^{-1} \bar{D}\_n, \\ D^k_n &=
(W^k_n)^{-1} \bar{D}\_{n-k+1} \\ D^{k-1}\_n, \quad \text{for \$k \geq
2\$}. \end{aligned} \$\$ We note that \\\bar{D}\_{n-k+1}\\ above denotes
the \\(n-k) \times (n-k+1)\\ version of the first difference matrix that
is defined in the second-to-last display.

The option `tf_weighting = TRUE` returns \\W^k_n D^k_n\\ where \\W^k_n\\
is the \\(n-k) \times (n-k)\\ diagonal matrix as described above. This
weighting is implicit in trend filtering, as explained in the help file
for
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md).
See also Section 6.1 of Tibshirani (2020) for further discussion.

**Note:** For multiplication of a given vector by \\D^k_n\\, instead of
forming \\D^k_n\\ with the current function and then carrying out the
multiplication, one should instead use
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md),
as this will be more efficient (both will be linear time, but the latter
saves the cost of forming any matrix in the first place).

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.1.

## See also

[`b_mat()`](https://glmgen.github.io/dspline/reference/b_mat.md) for
constructing the extended discrete derivative matrix, and
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md)
for multiplying by the discrete derivative matrix.

## Examples

``` r
d_mat(2, 1:10)
#> 8 x 10 sparse Matrix of class "dgCMatrix"
#>                                 
#> [1,] 1 -2  1  .  .  .  .  .  . .
#> [2,] .  1 -2  1  .  .  .  .  . .
#> [3,] .  .  1 -2  1  .  .  .  . .
#> [4,] .  .  .  1 -2  1  .  .  . .
#> [5,] .  .  .  .  1 -2  1  .  . .
#> [6,] .  .  .  .  .  1 -2  1  . .
#> [7,] .  .  .  .  .  .  1 -2  1 .
#> [8,] .  .  .  .  .  .  .  1 -2 1
d_mat(2, 1:10 / 10)
#> 8 x 10 sparse Matrix of class "dgCMatrix"
#>                                                     
#> [1,] 100 -200  100    .    .    .    .    .    .   .
#> [2,]   .  100 -200  100    .    .    .    .    .   .
#> [3,]   .    .  100 -200  100    .    .    .    .   .
#> [4,]   .    .    .  100 -200  100    .    .    .   .
#> [5,]   .    .    .    .  100 -200  100    .    .   .
#> [6,]   .    .    .    .    .  100 -200  100    .   .
#> [7,]   .    .    .    .    .    .  100 -200  100   .
#> [8,]   .    .    .    .    .    .    .  100 -200 100
d_mat(2, 1:10, row_idx = 2:5)
#> 4 x 10 sparse Matrix of class "dgCMatrix"
#>                             
#> [1,] . 1 -2  1  .  . . . . .
#> [2,] . .  1 -2  1  . . . . .
#> [3,] . .  .  1 -2  1 . . . .
#> [4,] . .  .  .  1 -2 1 . . .
```
