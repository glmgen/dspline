check_length <- function(x, n, op = c("==", ">=", "<=")) {
  args = as.character(as.list(match.call())[-1])
  op = match.arg(op)
  if (!(is.vector(x) && is.atomic(x))) {
    rlang::abort(sprintf("`%s` must be an atomic vector.", args[1]))
  }
  if (op == "==") {
    if (length(x) != n) {
      rlang::abort(sprintf("`length(%s)` must equal `%s`.", args[1], args[2]))
    }
  }
  else if (op == ">=") {
    if (length(x) < n) {
      rlang::abort(sprintf("`length(%s)` must be at least `%s`.", args[1], 
                           args[2]))
    }
  }
  else {
    if (length(x) > n) {
      rlang::abort(sprintf("`length(%s)` must be at most `%s`.", args[1], 
                           args[2]))
    }
  }
}

check_rows <- function(x, n) {
  args = as.character(as.list(match.call())[-1])
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    rlang::abort(sprintf("`%s` must be a matrix.", args[1]))
  }
  if (nrow(x) != n) {
    rlang::abort(sprintf("`nrow(%s)` must equal `%s`.", args[1], args[2]))
  }
}

check_cols <- function(x, n) {
  args = as.character(as.list(match.call())[-1])
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    rlang::abort(sprintf("`%s` must be a matrix.", args[1]))
  }
  if (ncol(x) != n) {
    rlang::abort(sprintf("`ncol(%s)` must equal `%s`.", args[1], args[2]))
  }
}

check_range <- function(x, rg) {
  args = as.character(as.list(match.call())[-1])
  if (!(min(x) >= min(rg) && max(x) <= max(rg))) {
    rlang::abort(sprintf("`%s` must lie in the range `%s.`", args[1], args[2]))
  }
}

check_sorted <- function(x) {
  args = as.character(as.list(match.call())[-1])
  if (!isTRUE(!is.unsorted(x))) {
    rlang::abort(sprintf("`%s` must be sorted in increasing order.", args[1]))
  }
}

check_nonneg_int <- function(x) {
  args = as.character(as.list(match.call())[-1])
  if (!(x >= 0 && round(x) == x)) {
    rlang::abort(sprintf("`%s` must be a nonnegative integer.", args[1]))
  }
}

check_no_nas <- function(x,
                         x_arg = rlang::caller_arg(x),
                         call = rlang::caller_call(),
                         ...) {
  if (anyNA(x)) {
    rlang::abort(sprintf("`%s` must not have any NAs.", x_arg), call = call,
                 ...)
  }
}
