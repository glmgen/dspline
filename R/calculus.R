#' Discrete differentiation
#'
#' Computes the discrete derivative of a function (or vector of function  
#' evaluations) of a given order, with respect to given design points, and
#' evaluated at a given query point.  
#'
#' @param f Function, or vector of function evaluations at `c(xd, x)`, the
#'   design points `xd` adjoined with the query point(s) `x`.
#' @param k Order for the discrete derivative calculation.
#' @param xd Design points. These must be sorted in increasing order.
#' @param x Query point(s).
#' @return Discrete derivative of `f` of order `k`, with respect to design
#'   points `xd`, evaluated at the query point(s) `x`.
#'
#' @details The discrete derivative of a function \eqn{f} of order \eqn{k}, with
#'   respect design points \eqn{x_1 < \ldots < x_n}, and evaluated at a query
#'   point \eqn{x}, is defined as
#'   \deqn{
#'   (\Delta^k_n f) (x) =
#'   \begin{cases}
#'   k! \cdot f[x_{i-k+1},\ldots,x_i,x] & \text{if $x \in (x_i,x_{i+1}]$, $i
#'   \geq k$} \\
#'   i! \cdot f[x_1,\ldots,x_i,x] & \text{if $x \in (x_i,x_{i+1}]$, $i < k$} \\
#'   f(x) & \text{if $x \leq x_1$},
#'   \end{cases}
#'   }
#'   where we take \eqn{x_{n+1} = \infty} for convenience. In other words, for
#'   "most" points \eqn{x > x_k}, we define \eqn{(\Delta^k_n f)(x)} in terms of
#'   a (scaled) divided difference of \eqn{f} of order \eqn{k}, where the
#'   centers are the \eqn{k} points immediately to the left of \eqn{x}, and
#'   \eqn{x} itself. Meanwhile, for "boundary" points \eqn{x \leq x_k}, we
#'   define \eqn{(\Delta^k_n f)(x)} to be a (scaled) divided difference of
#'   \eqn{f} of the highest possible order, where the centers are the points to
#'   the left of \eqn{x}, and \eqn{x} itself. For more discussion, including
#'   alternative representations for this discrete differentiation operator, see
#'   Section 3.1 of Tibshirani (2020).
#'
#' Lastly, an **important note:** for calculation of discrete derivatives at the
#'   design points themselves, which could be achieved by taking `x = xd` in the 
#'   current function, one instead should use `b_mat_mult()` or `d_mat_mult()`,
#'   as these will be more efficient.
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 3.1.
#' @export
discrete_deriv <- function(f, k, xd, x) {
  if (is.function(f)) f = sapply(c(xd, x), f)
  rcpp_discrete_deriv(f, k, xd, x)
}

#' Discrete integration
#'
#' Computes the discrete integral of a function (or vector of function
#' evaluations) of a given order, with respect to given design points, and 
#' evaluated at a given query point.
#'
#' @param f Function, or vector of function evaluations at `c(xd, x)`, the 
#'   design points `xd` adjoined with the query point(s) `x`.
#' @param k Order for the discrete integral calculation.
#' @param xd Design points. These must be sorted in increasing order.
#' @param x Query point(s).
#' @return Discrete integral of `f` of order `k`, with respect to design
#'   points `xd`, evaluated at the query point(s) `x`.
#'
#' @details The discrete integral of a function \eqn{f} of order \eqn{k}, with 
#'   respect design points \eqn{x_1 < \ldots < x_n}, and evaluated at a query
#'   point \eqn{x}, is defined as TODO
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and 
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 3.2.
#' @export
discrete_integ <- function(f, k, xd, x) {
  if (is.function(f)) f = sapply(c(xd, x), f)
  rcpp_discrete_integ(f, k, xd, x)
}
