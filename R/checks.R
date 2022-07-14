check_length <- function(x, n, op = c("==", ">=", "<=")) {
  args = as.character(as.list(match.call())[-1])
  op = match.arg(op)
  if (op == "==") {
    if (length(x) != n) {
      rlang::abort(sprintf("`length(%s)` must equal `%s`.", args[1], args[2]))
    }
  }
  else if (op == ">=") {
    if (length(x) < n) {
      rlang::abort(sprintf("`length(%s)` must be at least `%s`.", args[1], args[2]))
    }
  }
  else {
    if (length(x) > n) {
      rlang::abort(sprintf("`length(%s)` must be at most `%s`.", args[1], args[2]))
    }
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
  
