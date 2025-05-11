#' @rdname dot_functions
#' @name dot_functions
#' @title In-place computations
#'
#' @description Each "dot" function accepts arguments as in its "non-dot"
#'   counterpart, but peforms computations in place, overwriting the first input
#'   argument (which must be a vector of the appropriate length) with the
#'   desired output.
#'
#' @details These functions should not be used unless you are intentionally
#'   doing so for memory considerations and are nonetheless being very careful.
#'
#' An **important warning:** each "dot" function only works as expected if its
#'   first argument is passed in as a vector of numeric type. If the first
#'   argument is passed in as an integer vector, then since the output must (in
#'   general) be a numeric vector, it cannot be computed in place (Rcpp performs
#'   an implicit cast and copy when it converts this to NumericVector type for
#'   use in C++).
#'
#' Also, each "dot" function does not perform any error checking on its input
#'   arguments. Use with care. More details on the computations performed by
#'   individual functions are provided below.
#'
#' @section `.divided_diff()`:
#' Overwrites `f` with all lower-order divided differences: each element `f[i]`
#'   becomes the divided difference with respect to centers `z[1:i]`. See also
#'   [divided_diff()].
#'
#' @section `.b_mat_mult()`:
#' Overwrites `v` with `B %*% v`, where `B` is the extended discrete derivative
#'   matrix as returned by `b_mat()`. See also [b_mat_mult()].
#'
#' @section `.h_mat_mult()`:
#' Overwrites `v` with `H %*% v`, where `H` is the falling factorial basis
#'   matrix as returned by `h_mat()`. See also [h_mat_mult()].
#' @return None. These functions *overwrite* their input.
#' @examples
#' v = as.numeric(1:10) # Note: must be of numeric type
#' b_mat_mult(v, 1, 1:10) 
#' v
#' .b_mat_mult(v, 1, 1:10, FALSE, FALSE, FALSE)
#' v
NULL

#' @rdname dot_functions
#' @inheritParams divided_diff
#' @export
.divided_diff <- rcpp_dot_divided_diff

#' @rdname dot_functions
#' @inheritParams b_mat_mult
#' @export
.b_mat_mult <- rcpp_dot_b_mat_mult

#' @rdname dot_functions
#' @inheritParams h_mat_mult
#' @export
.h_mat_mult <- rcpp_dot_h_mat_mult
