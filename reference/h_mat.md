# Construct H matrix

Constructs the falling factorial basis matrix of a given order, with
respect to given design points.

## Usage

``` r
h_mat(k, xd, di_weighting = FALSE, col_idx = NULL)
```

## Arguments

- k:

  Order for the falling factorial basis matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- di_weighting:

  Should "discrete integration weighting" be used? Multiplication by
  such a weighted H gives discrete integrals at the design points; see
  details for more information. The default is `FALSE`.

- col_idx:

  Vector of indices, a subset of `1:n` where `n = length(xd)`, that
  indicates which columns of the constructed matrix should be returned.
  The default is `NULL`, which is taken to mean `1:n`.

## Value

Sparse matrix of dimension `length(xd)` by `length(col_idx)`.

## Details

The falling factorial basis matrix of order \\k\\, with respect to
design points \\x_1 \< \ldots \< x_n\\, is denoted \\H^k_n\\. It has
dimension \\n \times n\\, and its entries are defined as: \$\$
(H^k_n)\_{ij} = h^k_j(x_i), \$\$ where \\h^k_1, \ldots, h^k_n\\ are the
falling factorial basis functions, defined as: \$\$ \begin{aligned}
h^k_j(x) &= \frac{1}{(j-1)!} \prod\_{\ell=1}^{j-1}(x-x\_\ell), \quad
j=1,\ldots,k+1, \\ h^k_j(x) &= \frac{1}{k!} \prod\_{\ell=j-k}^{j-1}
(x-x\_\ell) \cdot 1\\x \> x\_{j-1}\\, \quad j=k+2,\ldots,n.
\end{aligned} \$\$

The matrix \\H^k_n\\ can also be constructed recursively, as follows. We
first define the \\n \times n\\ lower triangular matrix of 1s: \$\$ L_n
= \left\[\begin{array}{rrrr} 1 & 0 & \ldots & 0 \\ 1 & 1 & \ldots & 0 \\
\vdots & & & \\ 1 & 1 & \ldots & 1 \end{array}\right\], \$\$ and for \\k
\geq 1\\, define the \\n \times n\\ extended diagonal weight matrix
\\Z^k_n\\ to have first \\k\\ diagonal entries equal to 1 and last
\\n-k\\ diagonal entries equal to \\(x\_{i+k} - x_i) / k\\, \\i =
1,\ldots,n-k\\. The \\k\\th order falling factorial basis matrix is then
given by the recursion: \$\$ \begin{aligned} H^0_n &= L_n, \\ H^k_n &=
H^{k-1}\_n Z^k_n \left\[\begin{array}{cc} I_k & 0 \\ 0 & L\_{n-k}
\end{array}\right\], \quad \text{for \$k \geq 1\$}, \end{aligned} \$\$
where \\I_k\\ denotes the \\k \times k\\ identity matrix, and
\\L\_{n-k}\\ denotes the \\(n-k) \times (n-k)\\ lower triangular matrix
of 1s. For further details about this recursive representation, see
Sections 3.3 and 6.3 of Tibshirani (2020).

The option `di_weighting = TRUE` returns \\H^k_n Z^{k+1}\_n\\ where
\\Z^{k+1}\_n\\ is the \\n \times n\\ diagonal matrix as defined above.
This is connected to discrete integration as explained in the help file
for
[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md).
See also Section 3.3 of Tibshirani (2020) for more details.

Each basis function \\h^k_j\\, for \\j \geq k+2\\, has a single knot at
\\x\_{j-1}\\. The falling factorial basis thus spans \\k\\th degree
piecewise polynomials—discrete splines, in fact—with knots in
\\x\_{(k+1):(n-1)}\\. The dimension of this space is \\n-k-1\\ (number
of knots) \\+\\ \\k+1\\ (polynomial dimension) \\=\\ \\n\\. Setting the
argument `col_idx` appropriately allow one to form a basis matrix for a
discrete spline space corresponding to an arbitrary knot set \\T
\subseteq x\_{(k+1):(n-1)}\\. For more information, see Sections 4.1 and
8 of Tibshirani (2020).

**Note 1:** For computing the least squares projection onto a discrete
spline space defined by an arbitrary knot set \\T \subseteq
x\_{(k+1):(n-1)}\\, one should **not** use the falling factorial basis,
but instead use the discrete natural spline basis from
[`n_mat()`](https://glmgen.github.io/dspline/reference/n_mat.md), as the
latter has **much** better numerical properties in general. The help
file for
[`dspline_solve()`](https://glmgen.github.io/dspline/reference/dspline_solve.md)
gives more information.

**Note 2:** For multiplication of a given vector by \\H^k_n\\, one
should **not** form \\H^k_n\\ with the current function and then carry
out the multiplication, but instead use
[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md),
as the latter will be **much** more efficient (quadratic-time versus
linear-time).

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.3.

## See also

[`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md)
for multiplying by the falling factorial basis matrix and
[`h_eval()`](https://glmgen.github.io/dspline/reference/h_eval.md) for
constructing evaluations of the falling factorial basis at arbitrary
query points.

## Examples

``` r
h_mat(2, 1:10)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                               
#>  [1,] 1 .  .  .  .  .  . . . .
#>  [2,] 1 1  .  .  .  .  . . . .
#>  [3,] 1 2  1  .  .  .  . . . .
#>  [4,] 1 3  3  1  .  .  . . . .
#>  [5,] 1 4  6  3  1  .  . . . .
#>  [6,] 1 5 10  6  3  1  . . . .
#>  [7,] 1 6 15 10  6  3  1 . . .
#>  [8,] 1 7 21 15 10  6  3 1 . .
#>  [9,] 1 8 28 21 15 10  6 3 1 .
#> [10,] 1 9 36 28 21 15 10 6 3 1
h_mat(2, 1:10 / 10)
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
h_mat(2, 1:10, col_idx = 4:7)
#> 10 x 4 sparse Matrix of class "dgCMatrix"
#>                  
#>  [1,]  .  .  .  .
#>  [2,]  .  .  .  .
#>  [3,]  .  .  .  .
#>  [4,]  1  .  .  .
#>  [5,]  3  1  .  .
#>  [6,]  6  3  1  .
#>  [7,] 10  6  3  1
#>  [8,] 15 10  6  3
#>  [9,] 21 15 10  6
#> [10,] 28 21 15 10
```
