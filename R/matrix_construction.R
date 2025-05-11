#' Construct D matrix
#'
#' Constructs the discrete derivative matrix of a given order, with respect to
#' given design points.
#'
#' @param k Order for the discrete derivative matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param tf_weighting Should "trend filtering weighting" be used? This is a
#'   weighting of the discrete derivatives that is implicit in trend filtering;
#'   see details for more information. The default is `FALSE`.
#' @param row_idx Vector of indices, a subset of `1:(n-k)` where `n =
#'   length(xd)`, that indicates which rows of the constructed matrix should be
#'   returned. The default is `NULL`, which is taken to mean `1:(n-k)`.
#' @return Sparse matrix of dimension `length(row_idx)` by `length(xd)`.
#'
#' @details The discrete derivative matrix of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{D^k_n}. It has
#'   dimension \eqn{(n-k) \times n}, and is banded with bandwidth \eqn{k+1}. It
#'   can be constructed recursively, as follows. We first define the \eqn{(n-1)
#'   \times n} first difference matrix \eqn{\bar{D}_n}:
#'   \deqn{
#'   \bar{D}_n =
#'   \left[\begin{array}{rrrrrr}
#'   -1 & 1 & 0 & \ldots & 0 & 0 \\
#'   0 & -1 & 1 & \ldots & 0 & 0 \\
#'   \vdots & & & & & \\
#'   0 & 0 & 0 & \ldots & -1 & 1
#'   \end{array}\right],
#'   }
#'   and for \eqn{k \geq 1}, define the \eqn{(n-k) \times (n-k)} diagonal weight
#'   matrix \eqn{W^k_n} to have diagonal entries \eqn{(x_{i+k} - x_i) / k},
#'   \eqn{i = 1,\ldots,n-k}. The \eqn{k}th order discrete derivative matrix
#'   \eqn{D^k_n} is then given by the recursion:
#'   \deqn{
#'   \begin{aligned}
#'   D^1_n &= (W^1_n)^{-1} \bar{D}_n, \\
#'   D^k_n &= (W^k_n)^{-1} \bar{D}_{n-k+1} \, D^{k-1}_n,
#'   \quad \text{for $k \geq 2$}.
#'   \end{aligned}
#'   }
#'   We note that \eqn{\bar{D}_{n-k+1}} above denotes the \eqn{(n-k) \times
#'   (n-k+1)} version of the first difference matrix that is defined in the
#'   second-to-last display.
#'
#' The option `tf_weighting = TRUE` returns \eqn{W^k_n D^k_n} where \eqn{W^k_n}
#'   is the \eqn{(n-k) \times (n-k)} diagonal matrix as described above. This
#'   weighting is implicit in trend filtering, as explained in the help file for
#'   [d_mat_mult()]. See also Section 6.1 of Tibshirani (2020) for further
#'   discussion.
#'
#' **Note:** For multiplication of a given vector by \eqn{D^k_n}, instead of
#'   forming \eqn{D^k_n} with the current function and then carrying out the
#'   multiplication, one should instead use [d_mat_mult()], as this will be more
#'   efficient (both will be linear time, but the latter saves the cost of
#'   forming any matrix in the first place).
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.1.
#' @seealso [b_mat()] for constructing the extended discrete derivative matrix,
#'   and [d_mat_mult()] for multiplying by the discrete derivative matrix.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' d_mat(2, 1:10)
#' d_mat(2, 1:10 / 10)
#' d_mat(2, 1:10, row_idx = 2:5)
d_mat <- function(k, xd, tf_weighting = FALSE, row_idx = NULL) {
  check_nonneg_int(k)
  check_length(xd, k+1, ">=")
  check_sorted(xd)
  n = length(xd)
  if (is.null(row_idx)) row_idx = 1:(n-k)
  else check_range(row_idx, 1:(n-k))
  rcpp_b_mat(k, xd, tf_weighting, row_idx-1, TRUE)
}

#' Construct B matrix
#'
#' Constructs the extended discrete derivative matrix of a given order, with
#' respect to given design points.
#'
#' @param k Order for the extended discrete derivative matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param tf_weighting Should "trend filtering weighting" be used? This is a
#'   weighting of the discrete derivatives that is implicit in trend filtering;
#'   see details for more information. The default is `FALSE`.
#' @param row_idx Vector of indices, a subset of `1:n` where `n = length(xd)`,
#'   that indicates which rows of the constructed matrix should be returned. The
#'   default is `NULL`, which is taken to mean `1:n`.
#' @return Sparse matrix of dimension `length(row_idx)` by `length(xd)`.
#'
#' @details The extended discrete derivative matrix of order \eqn{k}, with
#'   respect to design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{B^k_n}.
#'   It has dimension \eqn{n \times n}, and is banded with bandwidth \eqn{k+1}.
#'   It can be constructed recursively, as follows. For \eqn{k \geq 1}, we first
#'   define the \eqn{n \times n} extended difference matrix \eqn{\bar{B}_{n,k}}:
#'   \deqn{
#'   \bar{B}_{n,k} =
#'   \left[\begin{array}{rrrrrrrrr}
#'   1 & 0 & \ldots & 0 & & & & \\
#'   0 & 1 & \ldots & 0 & & & & \\
#'   \vdots & & & & & & & \\
#'   0 & 0 & \ldots & 1 & & & & \\
#'   & & & -1 & 1 & 0 & \ldots & 0 & 0 \\
#'   & & & 0 & -1 & 1 & \ldots & 0 & 0 \\
#'   & & & \vdots & & & & & \\
#'   & & & 0 & 0 & 0 & \ldots & -1 & 1
#'   \end{array}\right]
#'   \begin{array}{ll}
#'   \left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
#'   \right\} & \hspace{-5pt} \text{$k$ rows} \\
#'   \left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
#'   \right\} & \hspace{-5pt} \text{$n-k$ rows}
#'   \end{array}.
#'   }
#'   We also define the \eqn{n \times n} extended diagonal weight matrix
#'   \eqn{Z^k_n} to have first \eqn{k} diagonal entries equal to 1 and last
#'   \eqn{n-k} diagonal entries equal to \eqn{(x_{i+k} - x_i) / k}, \eqn{i =
#'   1,\ldots,n-k}. The \eqn{k}th order extended discrete derivative matrix
#'   \eqn{B^k_n} is then given by the recursion:
#'   \deqn{
#'   \begin{aligned}
#'   B^1_n &= (Z^1_n)^{-1} \bar{B}_{n,1}, \\
#'   B^k_n &= (Z^k_n)^{-1} \bar{B}_{n,k} \, B^{k-1}_n,
#'   \quad \text{for $k \geq 2$}.
#'   \end{aligned}
#'   }
#'   We note that the discrete derivative matrix \eqn{D^k_n} from [d_mat()] is
#'   simply given by the last \eqn{n-k} rows of the extended matrix \eqn{B^k_n}.
#'
#' The option `tf_weighting = TRUE` returns \eqn{Z^k_n B^k_n} where \eqn{Z^k_n}
#'   is the \eqn{n \times n} diagonal matrix as described above. This weighting
#'   is implicit in trend filtering, as explained in the help file for
#'   [d_mat_mult()]. See also Sections 6.1 and 6.2 of Tibshirani (2020) for
#'   further discussion.
#'
#' **Note:** For multiplication of a given vector by \eqn{B^k_n}, instead of
#'   forming \eqn{B^k_n} with the current function and then carrying out the
#'   multiplication, one should instead use [b_mat_mult()], as this will be more
#'   efficient (both will be linear time, but the latter saves the cost of
#'   forming any matrix in the first place).
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.2.
#' @seealso [d_mat()] for constructing the discrete derivative matrix, and
#'   [b_mat_mult()] for multiplying by the extended discrete derivative matrix.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' b_mat(2, 1:10)
#' b_mat(2, 1:10 / 10)
#' b_mat(2, 1:10, row_idx = 4:7)
b_mat <- function(k, xd, tf_weighting = FALSE, row_idx = NULL) {
  check_nonneg_int(k)
  check_length(xd, k+1, ">=")
  check_sorted(xd)
  n = length(xd)
  if (is.null(row_idx)) row_idx = 1:n
  else check_range(row_idx, 1:n)
  rcpp_b_mat(k, xd, tf_weighting, row_idx-1, FALSE)
}

#' Construct H matrix
#'
#' Constructs the falling factorial basis matrix of a given order, with respect
#' to given design points.
#'
#' @param k Order for the falling factorial basis matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param di_weighting Should "discrete integration weighting" be used?
#'   Multiplication by such a weighted H gives discrete integrals at the design
#'   points; see details for more information. The default is `FALSE`.
#' @param col_idx Vector of indices, a subset of `1:n` where `n = length(xd)`,
#'   that indicates which columns of the constructed matrix should be
#'   returned. The default is `NULL`, which is taken to mean `1:n`.
#' @return Sparse matrix of dimension `length(xd)` by `length(col_idx)`.
#'
#' @details The falling factorial basis matrix of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{H^k_n}. It has
#'   dimension \eqn{n \times n}, and its entries are defined as:
#'   \deqn{
#'   (H^k_n)_{ij} = h^k_j(x_i),
#'   }
#'   where \eqn{h^k_1, \ldots, h^k_n} are the falling factorial basis functions,
#'   defined as:
#'   \deqn{
#'   \begin{aligned}
#'   h^k_j(x) &= \frac{1}{(j-1)!} \prod_{\ell=1}^{j-1}(x-x_\ell),
#'   \quad j=1,\ldots,k+1, \\
#'   h^k_j(x) &= \frac{1}{k!} \prod_{\ell=j-k}^{j-1} (x-x_\ell) \cdot
#'   1\{x > x_{j-1}\}, \quad j=k+2,\ldots,n.
#'   \end{aligned}
#'   }
#'
#' The matrix \eqn{H^k_n} can also be constructed recursively, as follows. We
#'   first define the \eqn{n \times n} lower triangular matrix of 1s:
#'   \deqn{
#'   L_n =
#'   \left[\begin{array}{rrrr}
#'   1 & 0 & \ldots & 0 \\
#'   1 & 1 & \ldots & 0 \\
#'   \vdots & & & \\
#'   1 & 1 & \ldots & 1
#'   \end{array}\right],
#'   }
#'   and for \eqn{k \geq 1}, define the \eqn{n \times n} extended diagonal
#'   weight matrix \eqn{Z^k_n} to have first \eqn{k} diagonal entries equal to 1
#'   and last \eqn{n-k} diagonal entries equal to \eqn{(x_{i+k} - x_i) / k},
#'   \eqn{i = 1,\ldots,n-k}. The \eqn{k}th order falling factorial basis matrix
#'   is then given by the recursion:
#'   \deqn{
#'   \begin{aligned}
#'   H^0_n &= L_n, \\
#'   H^k_n &= H^{k-1}_n Z^k_n
#'   \left[\begin{array}{cc}
#'   I_k & 0 \\
#'   0 & L_{n-k}
#'   \end{array}\right],
#'   \quad \text{for $k \geq 1$},
#'   \end{aligned}
#'   }
#'   where \eqn{I_k} denotes the \eqn{k \times k} identity matrix, and
#'   \eqn{L_{n-k}} denotes the \eqn{(n-k) \times (n-k)} lower triangular matrix
#'   of 1s.  For further details about this recursive representation, see
#'   Sections 3.3 and 6.3 of Tibshirani (2020).
#'
#' The option `di_weighting = TRUE` returns \eqn{H^k_n Z^{k+1}_n} where
#'   \eqn{Z^{k+1}_n} is the \eqn{n \times n} diagonal matrix as defined above.
#'   This is connected to discrete integration as explained in the help file for
#'   [h_mat_mult()]. See also Section 3.3 of Tibshirani (2020) for more details.
#'
#' Each basis function \eqn{h^k_j}, for \eqn{j \geq k+2}, has a single knot at
#'   \eqn{x_{j-1}}. The falling factorial basis thus spans \eqn{k}th degree
#'   piecewise polynomials---discrete splines, in fact---with knots in
#'   \eqn{x_{(k+1):(n-1)}}. The dimension of this space is \eqn{n-k-1} (number
#'   of knots) \eqn{+} \eqn{k+1} (polynomial dimension) \eqn{=} \eqn{n}. Setting
#'   the argument `col_idx` appropriately allow one to form a basis matrix for a
#'   discrete spline space corresponding to an arbitrary knot set \eqn{T
#'   \subseteq x_{(k+1):(n-1)}}. For more information, see Sections 4.1 and 8 of
#'   Tibshirani (2020).
#'
#' **Note 1:** For computing the least squares projection onto a discrete spline
#'   space defined by an arbitrary knot set \eqn{T \subseteq x_{(k+1):(n-1)}},
#'   one should **not** use the falling factorial basis, but instead use the
#'   discrete natural spline basis from [n_mat()], as the latter has **much**
#'   better numerical properties in general. The help file for [dspline_solve()]
#'   gives more information.
#'
#' **Note 2:** For multiplication of a given vector by \eqn{H^k_n}, one should
#'   **not** form \eqn{H^k_n} with the current function and then carry out the
#'   multiplication, but instead use [h_mat_mult()], as the latter will be
#'   **much** more efficient (quadratic-time versus linear-time).
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 6.3.
#' @seealso [h_mat_mult()] for multiplying by the falling factorial basis
#'   matrix and [h_eval()] for constructing evaluations of the falling factorial
#'   basis at arbitrary query points.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' h_mat(2, 1:10)
#' h_mat(2, 1:10 / 10)
#' h_mat(2, 1:10, col_idx = 4:7)
h_mat <- function(k, xd, di_weighting = FALSE, col_idx = NULL) {
  check_nonneg_int(k)
  check_length(xd, k+1, ">=")
  check_sorted(xd)
  n = length(xd)
  if (is.null(col_idx)) col_idx = 1:n
  else check_range(col_idx, 1:n)
  rcpp_h_mat(k, xd, di_weighting, col_idx-1)
}

#' Construct N matrix
#'
#' Constructs the discrete B-spline basis matrix of a given order, with respect
#' to given design points and given knot points.
#'
#' @param k Order for the discrete B-spline basis matrix. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param normalized Should the discrete B-spline basis vectors be normalized to
#'   attain a maximum value of 1 over the design points? The default is `TRUE`.
#' @param knot_idx Vector of indices, a subset of `(k+1):(n-1)` where `n =
#'   length(xd)`, that indicates which design points should be used as knot
#'   points for the discrete B-splines. Must be sorted in increasing order. The
#'   default is `NULL`, which is taken to mean `(k+1):(n-1)` as in the
#'   "canonical" discrete spline space (though in this case the returned N
#'   matrix will be trivial: it will be the identity matrix). See details.
#' @return Sparse matrix of dimension `length(xd)` by `length(knot_idx) + k +
#'   1`.
#'
#' @details The discrete B-spline basis matrix of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, and knot set \eqn{T \subseteq
#'   x_{(k+1):(n-1)}} is denoted \eqn{N^k_T}. It has dimension \eqn{(|T|+k+1)
#'   \times n}, and its entries are given by:
#'   \deqn{
#'   (N^k_T)_{ij} = \eta^k_j(x_i),
#'   }
#'   where \eqn{\eta^k_1, \ldots, \eta^k_m} are discrete B-spline (DB-spline)
#'   basis functions and \eqn{m = |T|+k+1}. As is suggested by their name, the
#'   DB-spline functions are linearly independent and span the space of
#'   discrete splines with knots at \eqn{T}. Each DB-spline \eqn{\eta^k_j} has a
#'   key local support property: it is supported on an interval containing at
#'   most \eqn{k+2} adjacent knots.
#'
#' The functions \eqn{\eta^k_1, \ldots, \eta^k_m} are, in general, not available
#'   in closed-form, and are defined by setting up and solving a sequence of
#'   locally-defined linear systems. For any knot set \eqn{T}, computation of
#'   the evaluations of all DB-splines at the design points can be done in
#'   \eqn{O(nk^2)} operations; see Sections 7, 8.2, and 8.3 of Tibshirani (2020)
#'   for details. The current function uses a sparse QR decomposition from the
#'   `Eigen::SparseQR` module in C++ in order to solve the local linear systems.
#'
#' When \eqn{T = x_{(k+1):(n-1)}}, the knot set corresponding to the "canonical"
#'   discrete spline space (spanned by the falling factorial basis functions
#'   \eqn{h^k_1, \ldots, h^k_n} whose evaluations make up \eqn{H^k_n}; see the
#'   help file for [h_mat()]), the DB-spline basis matrix, which we denote by
#'   \eqn{N^k_n}, is trivial: it equals the \eqn{n \times n} identity matrix,
#'   \eqn{N^k_n = I_n}. Therefore DB-splines are really only useful for knot
#'   sets \eqn{T} that are proper subsets of \eqn{x_{(k+1):(n-1)}}.
#'   Specification of the knot set \eqn{T} is done via the argument `knot_idx`.
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Sections 7, 8.2, and 8.3.
#' @seealso [h_eval()] for constructing evaluations of the discrete B-spline
#'   basis at arbitrary query points.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' n_mat(2, 1:10, knot_idx = c(3, 5, 7))
#' n_mat(2, 1:10, knot_idx = c(4, 6, 8))
n_mat <- function(k, xd, normalized = TRUE, knot_idx = NULL) {
  check_nonneg_int(k)
  check_length(xd, k+1, ">=")
  check_sorted(xd)
  n = length(xd)
  if (is.null(knot_idx)) {
    knot_idx = (k+1):(n-1)
  }
  else {
    check_range(knot_idx, (k+1):(n-1))
    check_sorted(knot_idx)
  }
  rcpp_n_mat(k, xd, normalized, knot_idx-1)
}
