# Multiply by H matrix

Multiplies a given vector by H, the falling factorial basis matrix of a
given order, with respect to given design points.

## Usage

``` r
h_mat_mult(v, k, xd, di_weighting = FALSE, transpose = FALSE, inverse = FALSE)
```

## Arguments

- v:

  Vector to be multiplied by H, the falling factorial basis matrix.

- k:

  Order for the falling factorial basis matrix. Must be \>= 0.

- xd:

  Design points. Must be sorted in increasing order, and have length at
  least `k+1`.

- di_weighting:

  Should "discrete integration weighting" be used? Multiplication by
  such a weighted H gives discrete integrals at the design points; see
  details for more information. The default is `FALSE`.

- transpose:

  Multiply by the transpose of H? The default is `FALSE`.

- inverse:

  Multiply by the inverse of H? The default is `FALSE`.

## Value

Product of falling factorial basis matrix H and the input vector `v`.

## Details

The falling factorial basis matrix of order \\k\\, with respect to
design points \\x_1 \< \ldots \< x_n\\, is denoted \\H^k_n\\. Its
entries are defined as: \$\$ (H^k_n)\_{ij} = h^k_j(x_i), \$\$ where
\\h^k_j\\ is the \\j\\th falling factorial basis function, as defined in
the help file for
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md). The
matrix \\H^k_n\\ can be constructed recursively as the product of
\\H^{k-1}\_n\\ and a diagonally-weighted cumulative sum matrix; see the
help file for
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md), or
Section 6.3 of Tibshirani (2020). Therefore, multiplication by \\H^k_n\\
or by its transpose can be performed in \\O(nk)\\ operations based on
iterated weighted cumulative sums. See Appendix D of Tibshirani (2020)
for details.

The option `di_weighting = TRUE` performs multiplication by \\H^k_n
Z^{k+1}\_n\\ where \\Z^{k+1}\_n\\ is an \\n \times n\\ diagonal matrix
whose first \\k+1\\ diagonal entries of \\Z^{k+1}\_n\\ are 1 and last
\\n-k-1\\ diagonal entries are \\(x\_{i+k+1} - x_i) / (k+1)\\, \\i =
1,\ldots,n-k-1\\. The connection to discrete integration is as follows:
multiplication of \\v = f(x\_{1:n})\\ by \\H^k_n Z^{k+1}\_n\\ gives
order \\k+1\\ discrete integrals (note the increment in order of
integration here) of \\f\\ at the points \\x\_{1:n}\\: \$\$ H^k_n
Z^{k+1}\_n v = (S^{k+1}\_n f)(x\_{1:n}). \$\$

Lastly, the matrix \\H^k_n\\ has a special **inverse relationship** to
the extended discrete derivative matrix \\B^{k+1}\_n\\ of degree
\\k+1\\; it satisfies: \$\$ H^k_n Z^{k+1}\_n B^{k+1}\_n = I_n, \$\$
where \\Z^{k+1}\_n\\ is the \\n \times n\\ diagonal matrix as described
above, and \\I_n\\ is the \\n \times n\\ identity matrix. This, combined
with the fact that the extended discrete derivative matrix has an
efficient recursive representation in terms of weighted differences,
means that multiplying by \\(H^k_n)^{-1}\\ or its transpose can be
performed in \\O(nk)\\ operations. See Section 6.2 and Appendix D of
Tibshirani (2020) for details.

## References

Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.2.

## See also

[`discrete_integ()`](https://glmgen.github.io/dspline/reference/discrete_integ.md)
for discrete integration at arbitrary query points, and
[`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md) for
constructing the falling factorial basis matrix.

## Examples

``` r
v = sort(runif(10))
as.vector(h_mat(2, 1:10) %*% v)
#>  [1]  0.01087246  0.07400980  0.23176991  0.74799459  1.92024386  4.06673045
#>  [7]  7.52670824 12.89401382 20.92940833 32.39414713
h_mat_mult(v, 2, 1:10) 
#>  [1]  0.01087246  0.07400980  0.23176991  0.74799459  1.92024386  4.06673045
#>  [7]  7.52670824 12.89401382 20.92940833 32.39414713
```
