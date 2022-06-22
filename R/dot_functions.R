#' @rdname dot_functions
#' @name dot_functions
#' @title In place computations
#' 
#' @description Each "dot" function accepts arguments as in its "non-dot"
#'   counterpart, but peforms computations in place, overwriting the first input
#'   argument (which must be a vector of the appropriate length) with the
#'   desired output.
#'
#' @details These functions should not be used unless you are intentionally
#'   doing so for memory considerations, and are very careful nonetheless.
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
#'   becomes the divided difference with respect to centers `z[1:i]`.
#'
#' @section `.d_mat_mult()`:
#'
#' @section `.b_mat_mult()`:
#'
#' @section `.h_mat_mult()`:
#' 
#'
NULL

#' @rdname dot_functions
#' @export
.divided_diff <- rcpp_dot_divided_diff
