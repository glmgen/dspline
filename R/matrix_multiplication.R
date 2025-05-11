#' Multiply by D matrix
#'
#' Multiplies a given vector by D, the discrete derivative matrix of a given
#' order, with respect to given design points.
#'
#' @param v Vector to be multiplied by D, the discrete derivative matrix.
#' @param k Order for the discrete derivative matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param tf_weighting Should "trend filtering weighting" be used? This is a
#'   weighting of the discrete derivatives that is implicit in trend filtering;
#'   see details for more information. The default is `FALSE`.
#' @param transpose Multiply by the transpose of D? The default is `FALSE`.
#' @return Product of the discrete derivative matrix D and the input vector `v`.
#'
#' @details The discrete derivative matrix of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{D^k_n}. It has
#'   dimension \eqn{(n-k) \times n}. Acting on a vector \eqn{v} of function
#'   evaluations at the design points, denoted \eqn{v = f(x_{1:n})}, it gives
#'   the discrete derivatives of \eqn{f} at the points \eqn{x_{(k+1):n}}:
#'   \deqn{
#'   D^k_n v = (\Delta^k_n f) (x_{(k+1):n}).
#'   }
#'   The matrix \eqn{D^k_n} can be constructed recursively as the product of a
#'   diagonally-weighted first difference matrix and \eqn{D^{k-1}_n}; see the
#'   help file for [d_mat()], or Section 6.1 of Tibshirani (2020). Therefore,
#'   multiplication by \eqn{D^k_n} or by its transpose can be performed in
#'   \eqn{O(nk)} operations based on iterated weighted differences. See Appendix
#'   D of Tibshirani (2020) for details.
#'
#' The option `tf_weighting = TRUE` performs multiplication by \eqn{W^k_n D^k_n}
#'   where \eqn{W^k_n} is a \eqn{(n-k) \times (n-k)} diagonal matrix with
#'   entries \eqn{(x_{i+k} - x_i) / k}, \eqn{i = 1,\ldots,n-k}. This weighting
#'   is implicit in trend filtering, as the penalty in the \eqn{k}th order trend
#'   filtering optimization problem (with optimization parameter \eqn{\theta})
#'   is \eqn{\|W^{k+1}_n D^{k+1}_n \theta\|_1}. Moreover, this is precisely the
#'   \eqn{k}th order total variation of the \eqn{k}th degree discrete spline
#'   interpolant \eqn{f} to \eqn{\theta}, with knots in \eqn{x_{(k+1):(n-1)}};
#'   that is, such an interpolant satisfies:
#'   \deqn{
#'   \mathrm{TV}(D^k f) = \|W^{k+1}_n D^{k+1}_n \theta\|_1,
#'   }
#'   where \eqn{D^k f} is the \eqn{k}th derivative of \eqn{f}. See Section
#'   9.1. of Tibshirani (2020) for more details.
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.1.
#' @seealso [discrete_deriv()] for discrete differentiation at arbitrary query
#'   points, [b_mat_mult()] for multiplying by the extended discrete derivative
#'   matrix, and [d_mat()] for constructing the discrete derivative matrix.
#' @export
#' @examples
#' v = sort(runif(10))
#' as.vector(d_mat(2, 1:10) %*% v)
#' d_mat_mult(v, 2, 1:10) 
d_mat_mult <- function(v, k, xd, tf_weighting = FALSE, transpose = FALSE) {
  check_nonneg_int(k)
  check_sorted(xd)
  check_length(xd, k+1, ">=")
  if (!transpose) check_length(v, length(xd))
  else check_length(v, length(xd) - k)
  rcpp_d_mat_mult(v, k, xd, tf_weighting, transpose)
}

#' Multiply by B matrix
#'
#' Multiplies a given vector by B, the extended discrete derivative matrix of a
#' given order, with respect to given design points.
#'
#' @param v Vector to be multiplied by B, the extended discrete derivative
#'   matrix.
#' @param k Order for the extended discrete derivative matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param tf_weighting Should "trend filtering weighting" be used? This is a
#'   weighting of the discrete derivatives that is implicit in trend filtering;
#'   see details for more information. The default is `FALSE`.
#' @param transpose Multiply by the transpose of B? The default is `FALSE`.
#' @param inverse Multiply by the inverse of B? The default is `FALSE`.
#' @return Product of the extended discrete derivative matrix B and the input
#'   vector `v`.
#'
#' @details The extended discrete derivative matrix of order \eqn{k}, with
#'   respect to design points \eqn{x_1 < \ldots < x_n}, is denoted
#'   \eqn{B^k_n}. It is square, having dimension \eqn{n \times n}. Acting on a
#'   vector \eqn{v} of function evaluations at the design points, denoted \eqn{v
#'   = f(x_{1:n})}, it gives the discrete derivatives of \eqn{f} at the points
#'   \eqn{x_{1:n}}:
#'   \deqn{
#'   B^k_n v = (\Delta^k_n f) (x_{1:n}).
#'   }
#'   The matrix \eqn{B^k_n} can be constructed recursively as the product of a
#'   diagonally-weighted first difference matrix and \eqn{B^{k-1}_n}; see the
#'   help file for [b_mat()], or Section 6.2 of Tibshirani (2020). Therefore,
#'   multiplication by \eqn{B^k_n} or by its transpose can be performed in
#'   \eqn{O(nk)} operations based on iterated weighted differences. See Appendix
#'   D of Tibshirani (2020) for details.
#'
#' The option `tf_weighting = TRUE` performs multiplication by \eqn{Z^k_n B^k_n}
#'   where \eqn{Z^k_n} is an \eqn{n \times n} diagonal matrix whose top left
#'   \eqn{k \times k} block equals the identity matrix and bottom right
#'   \eqn{(n-k) \times (n-k)} block equals \eqn{W^k_n}, the latter being a
#'   diagonal weight matrix that is implicit in trend filtering, as explained in
#'   the help file for [d_mat_mult()].
#'
#' Lastly, the matrix \eqn{B^k_n} has a special **inverse relationship** to the
#'   falling factorial basis matrix \eqn{H^{k-1}_n} of degree \eqn{k-1} with
#'   knots in \eqn{x_{k:(n-1)}}; it satisfies:
#'   \deqn{
#'   Z^k_n B^k_n H^{k-1}_n = I_n,
#'   }
#'   where \eqn{Z^k_n} is the \eqn{n \times n} diagonal matrix as described
#'   above, and \eqn{I_n} is the \eqn{n \times n} identity matrix. This,
#'   combined with the fact that the falling factorial basis matrix has an
#'   efficient recursive representation in terms of weighted cumulative sums,
#'   means that multiplying by \eqn{(B^k_n)^{-1}} or its transpose can be
#'   performed in \eqn{O(nk)} operations. See Section 6.3 and Appendix D of
#'   Tibshirani (2020) for details.
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.2.
#' @seealso [discrete_deriv()] for discrete differentiation at arbitrary query
#'   points, [d_mat_mult()] for multiplying by the discrete derivative matrix,
#'   and [b_mat()] for constructing the extended discrete derivative matrix.
#' @export
#' @examples
#' v = sort(runif(10))
#' as.vector(b_mat(2, 1:10) %*% v)
#' b_mat_mult(v, 2, 1:10) 
b_mat_mult <- function(v, k, xd, tf_weighting = FALSE, transpose = FALSE,
                       inverse = FALSE) {
  check_nonneg_int(k)
  check_sorted(xd)
  check_length(xd, k+1, ">=")
  check_length(v, length(xd))
  rcpp_b_mat_mult(v, k, xd, tf_weighting, transpose, inverse)
}

#' Multiply by H matrix
#'
#' Multiplies a given vector by H, the falling factorial basis matrix of a given
#' order, with respect to given design points.
#'
#' @param v Vector to be multiplied by H, the falling factorial basis matrix.
#' @param k Order for the falling factorial basis matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param di_weighting Should "discrete integration weighting" be used?
#'   Multiplication by such a weighted H gives discrete integrals at the
#'   design points; see details for more information. The default is `FALSE`.
#' @param transpose Multiply by the transpose of H? The default is `FALSE`.
#' @param inverse Multiply by the inverse of H? The default is `FALSE`.
#' @return Product of falling factorial basis matrix H and the input vector `v`.
#'
#' @details The falling factorial basis matrix of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{H^k_n}. Its entries
#'   are defined as:
#'   \deqn{
#'   (H^k_n)_{ij} = h^k_j(x_i),
#'   }
#'   where \eqn{h^k_j} is the \eqn{j}th falling factorial basis function, as
#'   defined in the help file for [h_mat()]. The matrix \eqn{H^k_n} can be
#'   constructed recursively as the product of \eqn{H^{k-1}_n} and a
#'   diagonally-weighted cumulative sum matrix; see the help file for [h_mat()],
#'   or Section 6.3 of Tibshirani (2020). Therefore, multiplication by
#'   \eqn{H^k_n} or by its transpose can be performed in \eqn{O(nk)} operations
#'   based on iterated weighted cumulative sums. See Appendix D of Tibshirani
#'   (2020) for details.
#'
#' The option `di_weighting = TRUE` performs multiplication by \eqn{H^k_n
#'   Z^{k+1}_n} where \eqn{Z^{k+1}_n} is an \eqn{n \times n} diagonal matrix
#'   whose first \eqn{k+1} diagonal entries of \eqn{Z^{k+1}_n} are 1 and last
#'   \eqn{n-k-1} diagonal entries are \eqn{(x_{i+k+1} - x_i) / (k+1)}, \eqn{i =
#'   1,\ldots,n-k-1}. The connection to discrete integration is as follows:
#'   multiplication of \eqn{v = f(x_{1:n})} by \eqn{H^k_n Z^{k+1}_n} gives order
#'   \eqn{k+1} discrete integrals (note the increment in order of integration
#'   here) of \eqn{f} at the points \eqn{x_{1:n}}:
#'   \deqn{
#'   H^k_n Z^{k+1}_n v = (S^{k+1}_n f)(x_{1:n}).
#'   }
#'
#' Lastly, the matrix \eqn{H^k_n} has a special **inverse relationship** to the
#'   extended discrete derivative matrix \eqn{B^{k+1}_n} of degree \eqn{k+1}; it
#'   satisfies:
#'   \deqn{
#'   H^k_n Z^{k+1}_n B^{k+1}_n = I_n,
#'   }
#'   where \eqn{Z^{k+1}_n} is the \eqn{n \times n} diagonal matrix as described
#'   above, and \eqn{I_n} is the \eqn{n \times n} identity matrix. This,
#'   combined with the fact that the extended discrete derivative matrix has an
#'   efficient recursive representation in terms of weighted differences, means
#'   that multiplying by \eqn{(H^k_n)^{-1}} or its transpose can be performed in
#'   \eqn{O(nk)} operations. See Section 6.2 and Appendix D of Tibshirani (2020)
#'   for details.
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.2.
#' @seealso [discrete_integ()] for discrete integration at arbitrary query
#'   points, and [h_mat()] for constructing the falling factorial basis matrix.
#' @export
#' @examples
#' v = sort(runif(10))
#' as.vector(h_mat(2, 1:10) %*% v)
#' h_mat_mult(v, 2, 1:10) 
h_mat_mult <- function(v, k, xd, di_weighting = FALSE, transpose = FALSE,
                       inverse = FALSE) {
  check_nonneg_int(k)
  check_sorted(xd)
  check_length(xd, k+1, ">=")
  check_length(v, length(xd))
  rcpp_h_mat_mult(v, k, xd, di_weighting, transpose, inverse)
}
