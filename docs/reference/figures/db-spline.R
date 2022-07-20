devtools::load_all()
n = 50
k = 2

set.seed(11)
#xd = sort(runif(n))
e = runif(n, -1/(2.5*(n+1)), 1/(2.5*(n+1)))
xd = 1:n/(n+1) + e
x = seq(0, 1, length = 5*n)
knot_idx = round(seq((k+1) + 5, (n-1) - 5, length = 4))
N1 = n_mat(k, xd, knot_idx = knot_idx)
N2 = n_eval(k, xd, x, knot_idx = knot_idx, N = N1)

png(file = "db-spline.png", width = 7, height = 3, units = "in", res = 500)
par(mar = rep(0.01, 4))
matplot(x, N2, type = "l", lty = 1, col = 1:ncol(N2), 
        xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
matplot(xd, N1, type = "p", pch = 19, col = 1:ncol(N1), add = TRUE)
abline(v = xd[knot_idx], lty = 2, lwd = 0.5, col = 8)
graphics.off()
