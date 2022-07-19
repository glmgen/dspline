devtools::load_all()
n = 100
k = 2

set.seed(11)
xd = sort(runif(n))
x = seq(0, 1, length=100)
knot_idx = round(seq((k+1) + 10, (n-1) - 10, length = 4))
N1 = n_mat(k, xd, knot_idx = knot_idx)
N2 = n_eval(k, xd, x, knot_idx = knot_idx, N = N1)

png(file = "db-spline.png", width = 7, height = 3, units = "in", res = 500)
par(mar = rep(0.01, 4))
matplot(x, N2, type = "l", lty = 1, col = 1:ncol(N2), 
        xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
matplot(xd, N1, type = "p", pch = 19, col = 1:ncol(N1), add = TRUE)
graphics.off()
