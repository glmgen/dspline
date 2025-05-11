#' Divided differencing
#'
#' Computes the divided difference of a function (or vector of function
#' evaluations) with respect to given centers.
#'
#' @param f Function, or vector of function evaluations at the centers.
#' @param z Centers for the divided difference calculation.
#' @return Divided difference of `f` with respect to centers `z`.
#'
#' @details The divided difference of a function \eqn{f} with respect to centers
#'   \eqn{z_1, \ldots, z_{k+1}} is defined recursively as:
#'   \deqn{
#'   f[z_1,\ldots,z_{k+1}] = \displaystyle
#'   \frac{f[z_2,\ldots,z_{k+1}] - f[z_1,\ldots,z_k]}{z_{k+1}-z_1},
#'   }
#'   with base case \eqn{f[z_1] = f(z_1)} (that is, divided differencing with
#'   respect to a single point reduces to function evaluation).
#'
#' A notable special case is when the centers are evenly-spaced, say, \eqn{z_i =
#'   z+ih}, \eqn{i=0,\ldots,k} for some spacing \eqn{h>0}, in which case the
#'   divided difference becomes a (scaled) forward difference, or equivalently a
#'   (scaled) backward difference,
#'   \deqn{
#'   k! \cdot f[z,\ldots,z+kh] = \displaystyle
#'   \frac{1}{h^k} (F^k_h f)(z) =
#'   \frac{1}{h^k} (B^k_h f)(z+kh),
#'   }
#'   where we use \eqn{F^k_h} and \eqn{B^k_v} to denote the forward and
#'   backward difference operators, respectively, of order \eqn{k} and with
#'   spacing \eqn{h}.
#'
#' @export
#' @examples
#' f = runif(4)
#' z = runif(4)
#' divided_diff(f[1], z[1])
#' f[1]
#' divided_diff(f[1:2], z[1:2])
#' (f[1]-f[2])/(z[1]-z[2])
#' divided_diff(f[1:3], z[1:3])
#' ((f[1]-f[2])/(z[1]-z[2]) - (f[2]-f[3])/(z[2]-z[3])) / (z[1]-z[3]) 
#' divided_diff(f, 1:4)
#' diff(f, diff = 3) / factorial(3)
divided_diff <- function(f, z) {
  if (is.function(f)) f = sapply(z, f)
  else check_length(f, length(z))
  rcpp_divided_diff(f, z)
}

#' Discrete differentiation
#'
#' Computes the discrete derivative of a function (or vector of function
#' evaluations) of a given order, with respect to given design points, and
#' evaluated at a given query point.
#'
#' @param f Function, or vector of function evaluations at `c(xd, x)`, the
#'   design points `xd` adjoined with the query point(s) `x`.
#' @param k Order for the discrete derivative calculation. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param x Query point(s).
#' @return Discrete derivative of `f` of order `k`, with respect to design
#'   points `xd`, evaluated at the query point(s) `x`.
#'
#' @details The discrete derivative operator of order \eqn{k}, with respect
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{\Delta^k_n}. Acting
#'   on a function \eqn{f}, and evaluated at a query point \eqn{x}, it yields:
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
#'   alternative representations for the discrete differentiation, see Section
#'   3.1 of Tibshirani (2020).
#'
#' **Note:** for calculating discrete derivatives at the design points
#'   themselves, which could be achieved by taking `x = xd` in the current
#'   function, one should instead use [b_mat_mult()] or [d_mat_mult()], as these
#'   will be more efficient (both will be linear-time, but the latter functions
#'   will be faster).
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 3.1.
#' @seealso [b_mat_mult()], [d_mat_mult()] for multiplication by the extended
#'   and non-extended discrete derivative matrices, giving discrete derivatives
#'   at design points.
#' @export
#' @examples
#' xd = 1:10 / 10
#' discrete_deriv(function(x) x^2, 1, xd, xd)
discrete_deriv <- function(f, k, xd, x) {
  check_nonneg_int(k)
  check_sorted(xd)
  check_length(xd, k+1, ">=")
  if (is.function(f)) f = sapply(c(xd, x), f)
  else check_length(f, length(xd) + length(x))
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
#' @param k Order for the discrete integral calculation. Must be >= 0.
#' @param xd Design points. Must be sorted in increasing order, and have length
#'   at least `k+1`.
#' @param x Query point(s).
#' @return Discrete integral of `f` of order `k`, with respect to design points
#'   `xd`, evaluated at the query point(s) `x`.
#'
#' @details The discrete integral operator of order \eqn{k}, with respect to
#'   design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{S^k_n}. It is the
#'   inverse operator to the discrete derivative operator \eqn{\Delta^k_n}, so
#'   that:
#'   \deqn{
#'   S^k_n \Delta^k_n = \Delta^k_n S^k_n = \mathrm{Id},
#'   } where
#'   \eqn{\mathrm{Id}} denotes the identity operator. It can also be represented
#'   in a more explicit form, as follows. Acting on a function \eqn{f} of order
#'   \eqn{k}, and evaluated at a query point \eqn{x}, it yields:
#'   \deqn{
#'   (S^k_n f)(x) =
#'   \begin{cases}
#'   \displaystyle
#'   \sum_{j=1}^k h^{k-1}_j(x) \cdot f(x_j) +
#'   \sum_{j=k+1}^i h^{k-1}_j(x) \cdot \frac{x_j-x_{j-k}}{k} \cdot f(x_j)
#'   + h^{k-1}_{i+1}(x) \cdot \frac{x-x_{i-k+1}}{k} \cdot f(x) & \\
#'   & \hspace{-75pt} \text{if $x \in (x_i,x_{i+1}]$, $i \geq k$} \\
#'   \displaystyle
#'   \sum_{j=1}^i h^{k-1}_j(x) \cdot f(x_j) \,+\, h^{k-1}_{i+1}(x) \cdot f(x)
#'   & \hspace{-75pt} \text{if $x \in (x_i,x_{i+1}]$, $i <  k$} \\
#'   f(x) & \hspace{-75pt} \text{if $x \leq x_1$},
#'   \end{cases}
#'   }
#'   where \eqn{h^{k-1}_1, \ldots, h^{k-1}_n} denote the falling factorial basis
#'   functions of degree \eqn{k-1}, with knots in \eqn{x_{k:(n-1)}}. The help
#'   file for [h_mat()] gives a definition of the falling factorial basis. It
#'   can be seen (due to the one-sided support of the falling factorial basis
#'   functions) that discrete integration at \eqn{x = x_i}, \eqn{i = 1,\ldots,n}
#'   is equivalent to multiplication by a weighted version of the falling
#'   factorial basis matrix. For more details, including an alternative
#'   recursive representation for discrete integration (that elucidates its
#'   relationship to discrete differentiation), see Section 3.2 of Tibshirani
#'   (2020).
#'
#' **Note:** for calculating discrete integrals at the design points themselves,
#'   which could be achieved by taking `x = xd` in the current function, one
#'   should instead use [h_mat_mult()] with `di_weighting = TRUE`, as this will
#'   be **much** more efficient (quadratic-time versus linear-time).
#'
#' @references Tibshirani (2020), "Divided differences, falling factorials, and
#'   discrete splines: Another look at trend filtering and related problems",
#'   Section 3.2.
#' @seealso [h_mat_mult()] for multiplication by the falling factorial basis
#'   matrix, giving a weighted analog of discrete integration at the design
#'   points.
#' @export
#' @examples
#' xd = 1:10 / 10
#' discrete_integ(function(x) 1, 1, xd, xd)
discrete_integ <- function(f, k, xd, x) {
  check_nonneg_int(k)
  check_sorted(xd)
  check_length(xd, k+1, ">=")
  if (is.function(f)) f = sapply(c(xd, x), f)
  else check_length(f, length(xd) + length(x))
  rcpp_discrete_integ(f, k, xd, x)
}
