#' Evaluate H basis
#'
#' Evaluates the falling factorial basis of a given order, with respect to given
#' design points, at arbitrary query points.
#'
#' @param k Order for the falling factorial basis. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param x Query points. Must be sorted in increasing order.
#' @param col_idx Vector of indices, a subset of `1:n` where `n = length(xd)`,
#'   that indicates which columns of the constructed matrix should be returned.
#'   The default is `NULL`, which is taken to mean `1:n`.
#' @return Sparse matrix of dimension `length(x)` by `length(col_idx)`.
#'
#' @details The falling factorial basis functions of order \eqn{k}, defined with
#'   respect to design points \eqn{x_1 < \ldots < x_n}, are denoted \eqn{h^k_1,
#'   \ldots, h^k_n}. For their precise definition and further references, see
#'   the help file for [h_mat()]. The current function produces a matrix of
#'   evaluations of the falling factorial basis at an arbitrary sequence of
#'   query points. For each query point \eqn{x}, this matrix has a corresponding
#'   row with entries:
#'   \deqn{
#'   h^k_j(x), \; j = 1, \ldots, n.
#'   }
#'
#' @seealso [h_mat()] for constructing evaluations of the falling factorial
#'   basis at the design points.
#' @export
#' @examples
#' xd = 1:10 / 10
#' x = 1:9 / 10 + 0.05
#' h_mat(2, xd)
#' h_eval(2, xd, x)
h_eval <- function(k, xd, x, col_idx = NULL) {
  check_nonneg_int(k)
  check_length(xd, k+1, ">=")
  check_sorted(xd)
  check_sorted(x)
  n = length(xd)
  if (is.null(col_idx)) col_idx = 1:n
  else check_range(col_idx, 1:n)
  rcpp_h_eval(k, xd, x, col_idx-1)
}

#' Evaluate N basis
#'
#' Evaluates the discrete B-spline basis of a given order, with respect to given
#' design points, evaluated at arbitrary query points.
#'
#' @param k Order for the discrete B-spline basis. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param x Query points. Must be sorted in increasing order.
#' @param normalized Should the discrete B-spline basis vectors be normalized to
#'   attain a maximum value of 1 over the design points? The default is `TRUE`.
#' @param knot_idx Vector of indices, a subset of `(k+1):(n-1)` where `n =
#'   length(xd)`, that indicates which design points should be used as knot
#'   points for the discrete B-splines. Must be sorted in increasing order. The
#'   default is `NULL`, which is taken to mean `(k+1):(n-1)`.
#' @param N Matrix of discrete B-spline evaluations at the design points. The
#'   default is `NULL`, which means that this is precomputed before constructing
#'   the matrix of discrete B-spline evaluations at the query points. If `N` is
#'   non-`NULL`, then the argument `normalized` will be ignored (as this would
#'   have only been used to construct N at the design points).
#' @return Sparse matrix of dimension `length(x)` by `length(knot_idx) + k + 1`.
#'
#' @details The discrete B-spline basis functions of order \eqn{k}, defined with
#'   respect to design points \eqn{x_1 < \ldots < x_n}, are denoted
#'   \eqn{\eta^k_1, \ldots, \eta^k_n}. For a discussion of their properties and
#'   further references, see the help file for [n_mat()]. The current function
#'   produces a matrix of evaluations of the discrete B-spline basis at an
#'   arbitrary sequence of query points. For each query point \eqn{x}, this
#'   matrix has a corresponding row with entries:
#'   \deqn{
#'   \eta^k_j(x), \; j = 1, \ldots, n.
#'   }
#'
#' Unlike the falling factorial basis, the discrete B-spline basis is not
#'   generally available in closed-form. Therefore, the current function (unlike
#'   [h_eval()]) will first check if it should precompute the evaluations of the
#'   discrete B-spline basis at the design points. If the argument `N` is
#'   non-`NULL`, then it will use this as the matrix of evaluations at the
#'   design points; if `N` is `NULL`, then it will call [n_mat()] to produce
#'   such a matrix, and will pass to this function the arguments `normalized`
#'   and `knot_idx` accordingly.
#'
#' After obtaining the matrix of discrete B-spline evaluations at the design
#'   points, the fast interpolation scheme from [dspline_interp()] is used to
#'   produce evaluations at the query points.
#'
#' @seealso [n_mat()] for constructing evaluations of the discrete B-spline
#'   basis at the design points.
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' xd = 1:10 / 10
#' x = 1:9 / 10 + 0.05
#' n_mat(2, xd, knot_idx = c(3, 5, 7))
#' n_eval(2, xd, x, knot_idx = c(3, 5, 7))
n_eval <- function(k, xd, x, normalized = TRUE, knot_idx = NULL, N = NULL) {
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
  check_sorted(x)
  if (is.null(N)) {
    rcpp_n_eval(k, xd, x, normalized, knot_idx-1)
  } else {
    if (ncol(N) != length(knot_idx)+k+1) {
      rlang::abort("`length(knot_idx) + k + 1` must equal `ncol(N)`.")
    }
    rcpp_n_eval_precomputed(k, xd, x, knot_idx-1, N)
  }
}

