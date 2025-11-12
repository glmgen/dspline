# Construct B matrix

Constructs the extended discrete derivative matrix of a given order,
with respect to given design points.

## Usage

``` r
b_mat(k, xd, tf_weighting = FALSE, row_idx = NULL)
```

## Arguments

- k:

  Order for the extended discrete derivative matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- tf_weighting:

  Should "trend filtering weighting" be used? This is a weighting of the
  discrete derivatives that is implicit in trend filtering; see details
  for more information. The default is `FALSE`.

- row_idx:

  Vector of indices, a subset of `1:n` where `n = length(xd)`, that
  indicates which rows of the constructed matrix should be returned. The
  default is `NULL`, which is taken to mean `1:n`.

## Value

Sparse matrix of dimension `length(row_idx)` by `length(xd)`.

## Details

The extended discrete derivative matrix of order \\k\\, with respect to
design points \\x_1 \< \ldots \< x_n\\, is denoted \\B^k_n\\. It has
dimension \\n \times n\\, and is banded with bandwidth \\k+1\\. It can
be constructed recursively, as follows. For \\k \geq 1\\, we first
define the \\n \times n\\ extended difference matrix \\\bar{B}\_{n,k}\\:
\$\$ \bar{B}\_{n,k} = \left\[\begin{array}{rrrrrrrrr} 1 & 0 & \ldots & 0
& & & & \\ 0 & 1 & \ldots & 0 & & & & \\ \vdots & & & & & & & \\ 0 & 0 &
\ldots & 1 & & & & \\ & & & -1 & 1 & 0 & \ldots & 0 & 0 \\ & & & 0 & -1
& 1 & \ldots & 0 & 0 \\ & & & \vdots & & & & & \\ & & & 0 & 0 & 0 &
\ldots & -1 & 1 \end{array}\right\] \begin{array}{ll}
\left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
\right\\ & \hspace{-5pt} \text{\$k\$ rows} \\
\left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
\right\\ & \hspace{-5pt} \text{\$n-k\$ rows} \end{array}. \$\$ We also
define the \\n \times n\\ extended diagonal weight matrix \\Z^k_n\\ to
have first \\k\\ diagonal entries equal to 1 and last \\n-k\\ diagonal
entries equal to \\(x\_{i+k} - x_i) / k\\, \\i = 1,\ldots,n-k\\. The
\\k\\th order extended discrete derivative matrix \\B^k_n\\ is then
given by the recursion: \$\$ \begin{aligned} B^1_n &= (Z^1_n)^{-1}
\bar{B}\_{n,1}, \\ B^k_n &= (Z^k_n)^{-1} \bar{B}\_{n,k} \\ B^{k-1}\_n,
\quad \text{for \$k \geq 2\$}. \end{aligned} \$\$ We note that the
discrete derivative matrix \\D^k_n\\ from
[`d_mat()`](https://glmgen.github.io/dspline/reference/d_mat.md) is
simply given by the last \\n-k\\ rows of the extended matrix \\B^k_n\\.

The option `tf_weighting = TRUE` returns \\Z^k_n B^k_n\\ where \\Z^k_n\\
is the \\n \times n\\ diagonal matrix as described above. This weighting
is implicit in trend filtering, as explained in the help file for
[`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md).
See also Sections 6.1 and 6.2 of Tibshirani (2020) for further
discussion.

**Note:** For multiplication of a given vector by \\B^k_n\\, instead of
forming \\B^k_n\\ with the current function and then carrying out the
multiplication, one should instead use
[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md),
as this will be more efficient (both will be linear time, but the latter
saves the cost of forming any matrix in the first place).

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.2.

## See also

[`d_mat()`](https://glmgen.github.io/dspline/reference/d_mat.md) for
constructing the discrete derivative matrix, and
[`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md)
for multiplying by the extended discrete derivative matrix.

## Examples

``` r
b_mat(2, 1:10)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                   
#>  [1,]  1  .  .  .  .  .  .  .  . .
#>  [2,] -1  1  .  .  .  .  .  .  . .
#>  [3,]  1 -2  1  .  .  .  .  .  . .
#>  [4,]  .  1 -2  1  .  .  .  .  . .
#>  [5,]  .  .  1 -2  1  .  .  .  . .
#>  [6,]  .  .  .  1 -2  1  .  .  . .
#>  [7,]  .  .  .  .  1 -2  1  .  . .
#>  [8,]  .  .  .  .  .  1 -2  1  . .
#>  [9,]  .  .  .  .  .  .  1 -2  1 .
#> [10,]  .  .  .  .  .  .  .  1 -2 1
b_mat(2, 1:10 / 10)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                      
#>  [1,]   1    .    .    .    .    .    .    .    .   .
#>  [2,] -10   10    .    .    .    .    .    .    .   .
#>  [3,] 100 -200  100    .    .    .    .    .    .   .
#>  [4,]   .  100 -200  100    .    .    .    .    .   .
#>  [5,]   .    .  100 -200  100    .    .    .    .   .
#>  [6,]   .    .    .  100 -200  100    .    .    .   .
#>  [7,]   .    .    .    .  100 -200  100    .    .   .
#>  [8,]   .    .    .    .    .  100 -200  100    .   .
#>  [9,]   .    .    .    .    .    .  100 -200  100   .
#> [10,]   .    .    .    .    .    .    .  100 -200 100
b_mat(2, 1:10, row_idx = 4:7)
#> 4 x 10 sparse Matrix of class "dgCMatrix"
#>                             
#> [1,] . 1 -2  1  .  . . . . .
#> [2,] . .  1 -2  1  . . . . .
#> [3,] . .  .  1 -2  1 . . . .
#> [4,] . .  .  .  1 -2 1 . . .
```
