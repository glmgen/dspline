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
#'   \eqn{z_1, \ldots, z_{k+1}} is defined recursively as
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
divided_diff <- function(f, z) {
  if (is.function(f)) f = sapply(z, f)
  rcpp_divided_diff(f, z)
}
