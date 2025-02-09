source("old-funs.R")
tol = 1e-8

test_that("Divided differences", {
  m = 4
  f = rnorm(m)
  z = rnorm(m)

  expect_equal(divided_diff(f[1], z[1]), f[1], tolerance = tol)
  expect_equal(divided_diff(f[1:2], z[1:2]), (f[1]-f[2])/(z[1]-z[2]), tolerance = tol)
  expect_equal(divided_diff(f[1:3], z[1:3]), ((f[1]-f[2])/(z[1]-z[2]) - (f[2]-f[3])/(z[2]-z[3]))/(z[1]-z[3]), tolerance = tol)
  expect_equal(divided_diff(f, 1:m), diff(f, diff = m-1) / factorial(m-1), tolerance = tol)  

  i = sample(m)
  expect_equal(divided_diff(f, z), divided_diff(f[i], z[i]), tolerance = tol)
})

test_that("In-place divided differences", {
  m = 4
  f = rnorm(m)
  z = rnorm(m)
  g = f + 1e-32
  .divided_diff(f, z)
  
  expect_equal(f[1], g[1], tolerance = tol)
  expect_equal(f[2], (g[1]-g[2])/(z[1]-z[2]), tolerance = tol)
  expect_equal(f[3], ((g[1]-g[2])/(z[1]-z[2]) - (g[2]-g[3])/(z[2]-z[3])) / (z[1]-z[3]), tolerance = tol)
  
  i = sample(m)
  expect_equal(f[m], divided_diff(g[i], z[i]), tolerance = tol)
})


test_that("Discrete derivatives", {
  n = 10
  k = 2
  f = function(x) sin(x)
  xd = sort(runif(n))
  x = c(0, 0.5, 0.7, 1.1)
  d1 = discrete_deriv(f, k, xd, x)

  d2 = c()
  for (x0 in x) {
    suppressWarnings({i = max(max(which(xd < x0)), 0)})
    kk = min(k, i)
    xx = c(xd[Seq(i-kk+1, i)], x0)
    if (kk == 0) d2 = c(d2, f(x0))
    else d2 = c(d2, getD(kk+1, kk, xx) %*% sapply(xx, f))
  }
  expect_equal(d1, d2, tolerance = tol)
  
  f1 = sapply(xd, f)
  d1 = discrete_deriv(f, k, xd, xd)
  d2 = as.numeric(getB(n, k, xd) %*% f1)
  expect_equal(d1, d2, tolerance = tol)
})

test_that("Discrete integrals", {
  n = 10
  k = 2
  f = function(x) sin(x)
  xd = sort(runif(n))
  x = c(xd, 1.1)
  f1 = sapply(x, f)
  Df = discrete_deriv(c(f1[1:n], f1), k, xd, x)
  f2 = discrete_integ(c(Df[1:n], Df), k, xd, x)
  expect_equal(f1, f2, tolerance = tol)
})

test_that("Construct D matrix", {
  n = 10
  k = 2
  xd = sort(runif(n))
  row_idx = sample(n-k, 4)
  D1 = d_mat(k, xd, row_idx = row_idx)
  D2 = getD(n, k, xd)[row_idx, ]
  expect_equal(D1, D2, tolerance = tol)

  D1 = d_mat(k, xd, tf_weight = TRUE, row_idx = row_idx)
  D2 = getDtil(n, k, xd)[row_idx, ]
  expect_equal(D1, D2, tolerance = tol)
})

test_that("Construct B matrix", {
  n = 10
  k = 2
  xd = sort(runif(n))
  row_idx = sample(n, 4)
  B1 = b_mat(k, xd, row_idx = row_idx)
  B2 = getB(n, k, xd)[row_idx, ]
  expect_equal(B1, B2, tolerance = tol)

  B1 = b_mat(k, xd, tf_weight = TRUE, row_idx = row_idx)
  B2 = getBtil(n, k, xd)[row_idx, ]
  expect_equal(B1, B2, tolerance = tol)
})

test_that("Construct H matrix", {
  n = 10
  k = 2
  col_idx = sample(n, 4)
  xd = sort(runif(n))
  H1 = as.matrix(h_mat(k, xd, col_idx = col_idx))
  H2 = getH(n, k, xd)[, col_idx]
  expect_equal(H1, H2, tolerance = tol, ignore_attr = TRUE)

  H = h_mat(k, xd, di_weight = TRUE)
  B = b_mat(k+1, xd)
  Id = as.matrix(H %*% B)
  expect_equal(Id, diag(n), tolerance = tol, ignore_attr = TRUE)
})


test_that("Multiply by D matrix", {
  n = 10
  k = 2
  xd = sort(runif(n))
  v = rnorm(n)
  a1 = as.numeric(d_mat(k, xd) %*% v)
  a2 = d_mat_mult(v, k, xd)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(Matrix::t(d_mat(k, xd)) %*% v[1:(n-k)])
  a2 = d_mat_mult(v[1:(n-k)], k, xd, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(d_mat(k, xd, tf_weight = TRUE) %*% v)
  a2 = d_mat_mult(v, k, xd, tf_weight = TRUE)
  expect_equal(a1, a2, tolerance = tol)
  
  a1 = as.numeric(Matrix::t(d_mat(k, xd, tf_weight = TRUE)) %*% v[1:(n-k)])
  a2 = d_mat_mult(v[1:(n-k)], k, xd, tf_weight = TRUE, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)
})


test_that("Multiply by B matrix", {
  n = 10
  k = 2
  xd = sort(runif(n))
  v = rnorm(n)
  a1 = as.numeric(b_mat(k, xd) %*% v)
  a2 = b_mat_mult(v, k, xd)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(Matrix::t(b_mat(k, xd)) %*% v)
  a2 = b_mat_mult(v, k, xd, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)
  
  a1 = solve(b_mat(k, xd), v)
  a2 = b_mat_mult(v, k, xd, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)
  
  a1 = solve(Matrix::t(b_mat(k, xd)), v)
  a2 = b_mat_mult(v, k, xd, transpose = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)
  
  a1 = as.numeric(b_mat(k, xd, tf_weight = TRUE) %*% v)
  a2 = b_mat_mult(v, k, xd, tf_weight = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(Matrix::t(b_mat(k, xd, tf_weight = TRUE)) %*% v)
  a2 = b_mat_mult(v, k, xd, tf_weight = TRUE, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(b_mat(k, xd, tf_weight = TRUE), v)
  a2 = b_mat_mult(v, k, xd, tf_weight = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(Matrix::t(b_mat(k, xd, tf_weight = TRUE)), v)
  a2 = b_mat_mult(v, k, xd, tf_weight = TRUE, transpose = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)
})

test_that("Multiply by H matrix", {
  n = 10
  k = 2
  xd = sort(runif(n))
  v = rnorm(n)
  a1 = as.numeric(h_mat(k, xd) %*% v)
  a2 = h_mat_mult(v, k, xd)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(Matrix::t(h_mat(k, xd)) %*% v)
  a2 = h_mat_mult(v, k, xd, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(h_mat(k, xd), v)
  a2 = h_mat_mult(v, k, xd, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(Matrix::t(h_mat(k, xd)), v)
  a2 = h_mat_mult(v, k, xd, transpose = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(h_mat(k, xd, di_weight = TRUE) %*% v)
  a2 = h_mat_mult(v, k, xd, di_weight = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = as.numeric(Matrix::t(h_mat(k, xd, di_weight = TRUE)) %*% v)
  a2 = h_mat_mult(v, k, xd, di_weight = TRUE, transpose = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(h_mat(k, xd, di_weight = TRUE), v)
  a2 = h_mat_mult(v, k, xd, di_weight = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)

  a1 = solve(Matrix::t(h_mat(k, xd, di_weight = TRUE)), v)
  a2 = h_mat_mult(v, k, xd, di_weight = TRUE, transpose = TRUE, inverse = TRUE)
  expect_equal(a1, a2, tolerance = tol)
})

test_that("Interpolation", {
  n = 10
  k = 2
  xd = sort(runif(n))
  x = seq(0, 1, length=100)
  v = xd^2 + 0.05 * rnorm(n)
  y1 = dspline_interp(v, k, xd, x, implicit = TRUE)
  y2 = dspline_interp(v, k, xd, x, implicit = FALSE)
  y3 = as.numeric(predict.tf.fitted(x, v, xd, k))
  expect_equal(y1, y2, tolerance = tol)
  expect_equal(y2, y3, tolerance = tol)
  
  # "Wide query" requesting that xd be echoed back as well:
  xd_wq = xd
  x_wq = c(x, xd)
  v_wq = v
  y1_wq = dspline_interp(v_wq, k, xd_wq, x_wq, implicit = TRUE)
  y2_wq = dspline_interp(v_wq, k, xd_wq, x_wq, implicit = FALSE)
  y3_wq = as.numeric(predict.tf.fitted(x_wq, v_wq, xd_wq, k))
  expect_equal(y1_wq, y2_wq, tolerance = tol)
  expect_equal(y2_wq, y3_wq, tolerance = tol)
  expect_equal(y3_wq, c(y3, v), tolerance = tol)
  
  # "Wide all" situation with design points also including both the "true" query
  # and "true" design points, with NAs in v:
  xd_wa = c(x, xd)
  x_wa = c(x, xd)
  v_wa = c(rep(NA, length(x)), v)
  o_wa = order(xd_wa)
  xd_wa <- xd_wa[o_wa]
  v_wa <- v_wa[o_wa]
  expect_error(dspline_interp(v_wa, k, xd_wa, x_wa),
               regexp = "`v` must not have any NAs.")
})

test_that("Newton interpolation", {
  k = 4
  n = k+1
  xd = sort(runif(n))
  v = xd^k
  x = seq(0, 1, length=100)
  y1 = dspline_interp(v, k, xd, x, implicit = TRUE)
  y2 = dspline_interp(v, k, xd, x, implicit = FALSE)
  y3 = newton_interp(v, xd, x)
  expect_equal(y1, y2, tolerance = tol)
  expect_equal(y2, y3, tolerance = tol)
})

test_that("Construct N matrix", {
  n = 50
  k = 2
  knot_idx = sort(sample((k+1):(n-1), 4))
  xd = sort(runif(n))
  N1 = n_mat(k, xd, normalized = FALSE, knot_idx = knot_idx)
  N2 = dbs.evals.sk(k, xd, knot_idx)
  expect_equal(N1, N2, tolerance = tol)

  H = h_mat(k, xd, col_idx = knot_idx+1)
  R1 = as.matrix(Matrix::qr.resid(Matrix::qr(N1), H))
  R2 = as.matrix(Matrix::qr.resid(Matrix::qr(N2), H))
  Ze = matrix(0, n, length(knot_idx))
  expect_equal(R1, Ze, tolerance = tol, ignore_attr = TRUE)
  expect_equal(R2, Ze, tolerance = tol, ignore_attr = TRUE)
})

test_that("Evaluate H basis", {
  n = 10
  k = 2
  col_idx = c(sample(1:(k+1), 2), sample((k+2):n, 4))
  xd = sort(runif(n))
  x = seq(0, 1, length=100)
  H = h_mat(k, xd, col_idx = col_idx)
  Hx = h_eval(k, xd, x, col_idx = col_idx)

  a = rnorm(length(col_idx))
  v = as.numeric(H %*% a)
  y1 = dspline_interp(v, k, xd, x)
  y2 = as.numeric(Hx %*% a)
  expect_equal(y1, y2, tolerance = tol)
})

test_that("Evaluate N basis", {
  n = 50
  k = 2
  knot_idx = sort(sample((k+1):(n-1), 4))
  xd = sort(runif(n))
  x = seq(0, 1, length=100)

  # n_mat already tested against dbs.evals.sk previously
  N_xd = n_mat(k, xd, knot_idx = knot_idx)
  N_x0 = as.matrix(n_eval(k, xd, x, knot_idx = knot_idx))
  N_x1 = as.matrix(n_eval(k, xd, x, knot_idx = knot_idx, N = N_xd))
  N_x2 = matrix(0, length(x), ncol(N_x1))
  for (j in 1:ncol(N_x2)) N_x2[,j] = dspline_interp(N_xd[,j], k, xd, x)
  expect_equal(N_x0, N_x1, tolerance = tol, ignore_attr = TRUE)
  expect_equal(N_x0, N_x2, tolerance = tol, ignore_attr = TRUE)
})

test_that("Projection", {
  n = 50
  k = 2
  knot_idx = sort(sample((k+1):(n-1), 20))
  xd = sort(runif(n))
  v = rnorm(n)
  
  obj1 = dspline_solve(v, k, xd, knot_idx = knot_idx, basis = "N")
  obj2 = dspline_solve(v, k, xd, knot_idx = knot_idx, basis = "B")
  obj3 = dspline_solve(v, k, xd, knot_idx = knot_idx, basis = "H")

  expect_equal(obj2$sol, obj3$sol, tolerance = tol)
  expect_equal(sum((v - obj1$fit)^2), sum((v - obj2$fit)^2), tolerance = tol)
  expect_equal(sum((v - obj1$fit)^2), sum((v - obj3$fit)^2), tolerance = tol)
})
