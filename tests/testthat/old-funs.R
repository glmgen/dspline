Seq = function(a,b) {
  if (a<=b) return(a:b)
  else return(integer(0))
}

Sign = function(x) {
  return(ifelse(x==0, 1, sign(x)))
}

Curve = function(f, a=0, b=1, N=100, ...) {
  x = seq(a,b,length=N)
  y = numeric(N)
  for (i in 1:N) y[i] = f(x[i])
  plot(x, y, type="l", ...)
}

symmetrize = function(mat, lower=TRUE) {
  tmat = t(mat)
  if (lower) tmat[lower.tri(tmat)] = mat[lower.tri(mat)]
  else tmat[upper.tri(tmat)] = mat[upper.tri(mat)]
  return(tmat)
}

cond.num = function(A) {
  s = eigen(A,symmetric=TRUE,only.values=TRUE)$val
  return(max(s)/min(s))
}

clockwise90 = function(a) { 
  return(t(a[nrow(a):1,]))
}

enlist = function(...) {
  result = list(...)
  n = result[[1]]
  if ((nargs() == 1) & is.character(n)) {
    result = as.list(seq(n))
    names(result) = n
    for (i in n) result[[i]] = get(i)
  }
  else {
    n = sys.call()
    n = as.character(n)[-1]
    n2 = names(result)
    if (!is.null(n2)) {
      which = n2 != ""
      n[which] = n2[which]
    }
    names(result) = n
  }
  return(result)
}

## Falling factorial functions

getT = function(n, k, x=1:n/n) {
  T = x
  if (k==0) T = T[-1]
  else if (k%%2==1) T = T[-c(1:((k+1)/2),(n-(k-1)/2):n)]  
  else if (k%%2==0) T = T[-c(1:((k+2)/2),(n-(k-2)/2):n)]
  return(T)
}

getG = function(n, k, x=1:n/n) {
  T = getT(n,k,x)
  G = matrix(0,n,n)
  G[,1] = rep(1,n)
  for (j in 2:n) {
    if (j<=k+1) {
      G[,j] = x^(j-1)/factorial(j-1)
    }
    else {
      G[,j] = pmax(x-T[j-k-1],0)^k/factorial(k)
    }
  }
  return(G)
}

getH = function(n, k, x=1:n/n) {
  return(evalH(x,k,(k+1):(n-1),x))
}

# Note: jj are the indices of the active knots in the x vector
evalH = function(x0, k, jj, x=1:n/n) {
  n0 = length(x0)
  m = length(jj)
  N = m+k+1
  H = matrix(0,n0,N)
  H[,1] = rep(1,n0)
  
  for (j in Seq(2,k+1)) {
    for (i in 1:n0) {
      H[i,j] = prod(x0[i]-x[1:(j-1)])/factorial(j-1)
    }
  }
  for (l in 1:length(jj)) {
    j = jj[l] 
    for (i in 1:n0) {
      if (x0[i] > x[j]) {
        if (k==0) H[i,l+k+1] = 1
        else H[i,l+k+1] = prod(x0[i]-x[(j-k+1):j])/factorial(k)
      }
    }
  }
  return(H)
}

library(Matrix)

Id = function(n) {
  return(bandSparse(n, k=0, diagonals=list(rep(1,n))))
}

# These are the "right" ones based on divided diffs
getB = function(n, k, x=1:n/n) {
  I = Id(n)
  D = bandSparse(n, k=c(-1,0), diagonals=list(rep(-1,n-1), rep(1,n)))
  B = I
  for (i in Seq(1,k)) {
    wts = c(rep(1,i), i/(x[(i+1):n]-x[1:(n-i)]))
    M = bdiag(I[Seq(1,i-1),Seq(1,i-1)], D[1:(n-i+1),1:(n-i+1)])
    B = (wts * M) %*% B
  }
  return(B)
}

getD = function(n, k, x=1:n/n) {
  getB(n,k,x)[-(1:k),]
}

# These were the ones used in previous TF papers
getBtil = function(n, k, x=1:n/n) {
  c(rep(1,k),(x[(k+1):n]-x[1:(n-k)])/k) * getB(n,k,x)
}

getDtil = function(n, k, x=1:n/n) {
  getBtil(n,k,x)[-(1:k),]
}

## Prediction functions

predict.tf.fitted = function(newx,beta,x,k) {
  n = length(x)
  B = getBtil(n,k+1,x)
  alpha = B %*% beta
  return(predict.tf.basis(newx,alpha,x,k))
}

predict.tf.basis = function(newx,alpha,x,k) {
  n = length(x)
  jj = (k+1):(n-1)
  H = evalH(newx,k,jj,x)
  return(H %*% alpha)
}

predict.tf.basis.sk = function(newx,alpha,x,k,jj) {
  H = evalH(newx,k,jj,x)
  return(H %*% alpha)
}

predict.tf.deriv = function(newx,beta,x,k,tol=1e-6) {
  i = which(abs(newx-x) < tol)
  if (length(i) > 0) return(median(beta[i]))

  if (max(x) <= newx) j0 = n
  else j0 = min(which(x > newx))

  if (j0 < k+1) {
    x0 = c(x[1:j0],newx,x[(j0+1):(k+1)])
    D0 = as.numeric(getD(k+2,k+1,x0))
    return(-sum(D0[-(j0+1)]*beta[1:(k+1)]) / D0[j0+1])
  }
  else {
    x0 = c(x[(j0-k):j0],newx)
    D0 = as.numeric(getD(k+2,k+1,x0))
    return(-sum(D0[-(k+2)]*beta[(j0-k):j0]) / D0[k+2])
  }
}

## B-spline functions

tpf = function(z, k) {
  return(function(x) (x-z)^k * (x>z))
}

# Assumes length(z) >= 1 (so that k >= 1)
tff = function(z) {
  return(Vectorize(function(x) prod(x-z) * (x>min(z))))
}

tnp = function(z) {
  return(Vectorize(function(x) prod(x-z) * (x>max(z))))
}

bs0 = function(z) {
  return(function(x) ifelse(x > z[1] && x <= x[2]))
}

bs = function(z) {
  k = length(z) - 2
  D = getD(k+2,k+1,z)
  funs = sapply(z, function(z0,k) return(tpf(z0,k)), k)
  return(Vectorize(function(x0) 
    return((-1)^(k+1) * (z[k+2]-z[1]) * (1/factorial(k+1)) *
             as.numeric(D %*% sapply(funs, function(f) f(x0))))))
}

dbs.v = function(z, v) {
  k = length(z) - 2
  if (k==0) return(bs(z)) # Just avoid indexing issues
  D = getD(k+2,k+1,z)
  funs = sapply(z, function(z0,k,v) return(tff(z0 + 0:(k-1)*v)), k, v)
  return(Vectorize(function(x0) 
    return((-1)^(k+1) * (z[k+2]-z[1]) * (1/factorial(k+1)) *
           as.numeric(D %*% sapply(funs, function(f) f(x0))))))
}

dbs = function(z, x) {
  k = length(z) - 2
  if (k==0) return(bs(z)) # Just avoid indexing issues
  D = getD(k+2,k+1,z) 
  v = rep(0,n)
  for (i in 1:n) {
    if (x[i] <= z[1]) v[i] = 0
    else if (x[i] >= z[k+2]) v[i] = 0
    else {
      vals = tnp(x[(i-k+1):i])(z)
      v[i] = (z[k+2]-z[1]) * (1/factorial(k+1)) *
        as.numeric(D %*% vals)
    }
  }
  return(dbs.fun(k,x,v))
}

dbs.fun = function(k, x, v) {
  return(Vectorize(function(newx) predict.tf.deriv(newx,v,x,k)))
}

dbs.fun.sk = function(k, x, alpha, jj) {
  return(Vectorize(function(newx) predict.tf.basis.sk(newx,alpha,x,k,jj)))
}

dbs.Amat.sk = function(k, x, jj, xright=NULL) {
  n = length(x)
  m = length(jj)
  N = m+k+1
  if (is.null(xright)) {
    dx = max(diff(x))
    xright = max(x) + (1:(k+2))*dx
  }
  J = c(jj, n:(n+k))
  X = c(x, xright)
  
  Amat = matrix(0,N+k+1,N) 
  for (i in 1:(k+1)) {
    ii = c(Seq(1,i), Seq(J[i]-k+1,J[i]+1))
    H = evalH(X[ii],k,J[1:i],X)
    y = rep(0,i+k+1); y[i] = 1
    alpha = lsfit(H,y,intercept=FALSE)$coef
    Amat[1:(k+1),i] = alpha[1:(k+1)]      # poly coefs
    Amat[k+1+(1:i),i] = alpha[-(1:(k+1))] # pp coefs
  }
  
  for (i in 1:m) {
    ii = c(Seq(J[i]-k, J[i]), J[i+1],
           Seq(J[i+k+1]-k+1, J[i+k+1]+1))
    H = evalH(X[ii],k,J[i:(i+k+1)],X)
    y = rep(0,2*k+3); y[k+2] = 1
    alpha = lsfit(H,y,intercept=FALSE)$coef
    Amat[1:(k+1),i+k+1] = alpha[1:(k+1)]            # poly coefs
    Amat[k+1+(i:(i+k+1)),i+k+1] = alpha[-(1:(k+1))] # pp coefs
  }
  
  return(Amat[1:N,]) # cut off boundary knots      
}

dbs.Bmat = function(k, x, a=NULL, b=NULL, xleft=NULL) {
  n = length(x)
  if (is.null(a) || is.null(b) || is.null(xleft)) {
    dx = min(diff(x))
    a = min(x) - dx
    b = max(x) + dx
    xleft = a + Seq(-k+1,-1)*dx
  }
  X = Vectorize(function(i) {
    if (1 <= i && i <= n) return(x[i])
    if (i==0) return(a)
    if (i==n+1) return(b)
    if (-k+1 <= i && i <= -1) return(xleft[i+k])
    else return(NA)
  })
  bvec = rep(1,n)
  if (k >= 1) {
    for (i in 1:n) {
      v = X(Seq(i-k+1,i))
      z = c(X(i-k), v, X(i+1))
      bvec[i] = (z[k+2]-z[1]) * (1/factorial(k+1)) *
        getD(k+2,k+1,z)[k+2] * tnp(v)(X(i+1))
    }
  }
  return(bandSparse(n, k=0, diagonals=list(bvec)))
}

# Faster and more stable way of evaluating DB-splines at the design
# points (in the sparse knot case)

dbs.evals.sk = function(k, x, jj, xright=NULL, tol=1e-10) {
  n = length(x)
  m = length(jj)
  N = m+k+1
  if (is.null(xright)) {
    dx = max(diff(x))
    xright = max(x) + (1:(k+1))*dx
  }
  J = c(jj, n:(n+k))
  X = c(x, xright)

  evals = Matrix(0,n+k+1,N+k+1)
  for (i in 1:(k+1)) {
    I = 1:(J[i]+1)
    jj = J[1:i]  
    z = rep(0,length(I)); z[i] = 1
    evals[I,i] = dbs.local.interp.left(k,jj,i,X[I],z,tol)
  }
  
  for (i in 1:m) {
    I = Seq(J[i]-k,J[i+k+1]+1)
    jj = J[i:(i+k+1)]-J[i]+k+1
    z = rep(0,length(I)); z[jj[2]] = 1
    evals[I,i+k+1] = dbs.local.interp(k,jj,X[I],z,tol)
  }
  
  return(evals[1:n,1:N]) # cut off boundary points and knots
}

# Local interpolation function for "left" DB-splines (first k+1)
# Note: jj are the indices of the active knots in the x vector
# Also: y is the observation vector, and we want to interpolate
# between its ith element and last knot

dbs.local.interp.left = function(k, jj, i, x, y, tol=1e-10) {
  n = length(x)
  B = getB(n,k+1,x)

  # Interpolate between point i+1 and last knot
  if (jj[i]-k+1 > i+1) {
    pts.ind = Seq(i+1,jj[i]-k)
    row.ind = setdiff((k+2):jj[i],jj[Seq(1,i-1)]+1)
    col.ind = 1:max(row.ind)
    I1 = pts.ind; I2 = setdiff(col.ind,I1)
    B1 = suppressWarnings(Matrix(B[row.ind,I1],ncol=length(I1)))
    B2 = suppressWarnings(Matrix(B[row.ind,I2],ncol=length(I2)))
    y[I1] = -solve(qr(B1,tol=tol), as.numeric(B2 %*% y[I2]))
  }

  return(y)
}

# Local interpolation function for DB-splines (last r functions)
# Note: jj are the indices of the active knots in the x vector
# Also: y is the observation vector, and we want to interpolate
# between first and second knot, and second knot and last knot

dbs.local.interp = function(k, jj, x, y, tol=1e-10) {
  n = length(x)
  B = getB(n,k+1,x)
  
  # Interpolate between first knot and second knot
  if (jj[2] > jj[1]+1) {
    # Grab a local part of discrete deriv matrix
    pts.ind = Seq(jj[1]+1,jj[2]-1)
    row.ind = pts.ind + 1
    col.ind = Seq(min(row.ind)-k-1, max(row.ind))
    # Solve for unknown function evaluations 
    I1 = pts.ind; I2 = setdiff(col.ind,I1)
    B1 = suppressWarnings(Matrix(B[row.ind,I1],ncol=length(I1)))
    B2 = suppressWarnings(Matrix(B[row.ind,I2],ncol=length(I2)))
    y[I1] = -solve(qr(B1,tol=tol), as.numeric(B2 %*% y[I2]))
  }
  
  # Interpolate between second knot and last knot
  if (jj[k+2]-k+1 > jj[2]+1) {
    # Grab a local part of discrete deriv matrix
    pts.ind = Seq(jj[2]+1,jj[k+2]-k)
    knt.ind = jj[jj %in% pts.ind]
    reg.ind = setdiff(pts.ind,knt.ind)
    # Row indices are tricky to compute. Insight: for each knot that's in 
    # jj[2]+1 to jj[k+2]-k, there's a point in jj[k+2]-k+1 to jj[k+2]-1
    # that's a regular point (not a knot); this is because the total number
    # of knots in jj[2]+1 to jj[k+2]-1 must be exactly k+1
    ii = Seq(jj[k+2]-k+1,jj[k+2]-1)
    row.ind = c(reg.ind, ii[!(ii %in% jj)]) + 1
    col.ind = Seq(min(row.ind)-k-1, max(row.ind))
    # Solve for unknown function evaluations 
    I1 = pts.ind; I2 = setdiff(col.ind,I1)
    B1 = suppressWarnings(Matrix(B[row.ind,I1],ncol=length(I1)))
    B2 = suppressWarnings(Matrix(B[row.ind,I2],ncol=length(I2)))
    y[I1] = -solve(qr(B1,tol=tol), as.numeric(B2 %*% y[I2]))
  }
  
  return(y)
}

## Discrete spline smoothing

ff_deriv = function(z, d, knot=TRUE) {
  k = length(z)
  deriv_str = paste0(paste(paste0("(x-", letters[1:k]),collapse=")*"),")")
  deriv_expr = parse(text=deriv_str)
  for (i in Seq(1,d)) {
    deriv_expr = D(deriv_expr, "x")
  }
  list_str = paste0(paste0("z[",1:k),"]")
  list_str = paste(paste(letters[1:k], list_str, sep="="), collapse=", ")
  list_str = paste0("list(x=x, ", list_str, ")")
  if (!knot) {
    return(function(x) {
      eval(deriv_expr, envir=eval(parse(text=list_str))) / factorial(k)
    })
  }
  else {
    return(function(x) {
      eval(deriv_expr, envir=eval(parse(text=list_str))) / factorial(k) *
        (x > max(z))
    })
  }
}

getM = function(n, m, x=1:n/n, a=x[1], b=x[n], del=1e-8) {
  k = 2*m-1
  M = matrix(0,n,n)
  
  for (i in (m+1):n) {
    if (i <= k+1) zi = x[1:(i-1)]
    else zi = x[(i-k):(i-1)]
    for (j in (m+1):i) {
      if (j <= k+1) zj = x[1:(j-1)]
      else zj = x[(j-k):(j-1)]
        
      # First min(m,i) - 1 terms
      for (l in Seq(1,min(i-m,m)-1)) {
        hi_deriv = ff_deriv(zi,m+l-1,knot=i>k+1)
        hj_deriv = ff_deriv(zj,m-l,knot=j>k+1)
  
        M[i,j] = M[i,j] + (-1)^(l-1) * (
          hi_deriv(b)*hj_deriv(b) -
          hi_deriv(x[i-1]+del)*hj_deriv(x[i-1]+del) +
          hi_deriv(x[i-1]-del)*hj_deriv(x[i-1]-del) -
          hi_deriv(a)*hj_deriv(a))
      }
      
      # Last term
      if (i <= k+1) {
        hj_deriv = ff_deriv(zj,2*m-i,knot=j>k+1)
        M[i,j] = M[i,j] + (-1)^(i-m-1) * (
          hj_deriv(b) - hj_deriv(a))
      }
      else {
        hj_deriv = ff_deriv(zj,0,knot=j>k+1)
        M[i,j] = M[i,j] + (-1)^(m-1) * 
          (hj_deriv(b) - hj_deriv(x[i-1]+del))
      }
    }
  }

  return(symmetrize(M))
}

getV = function(n, m, x=1:n/n, a=x[1], b=x[n]) {
  M = getM(n,m,x,a,b)
  Btil = getBtil(n,2*m,x)
  Htil = getH(n,m-1,x) %*% diag(c(rep(1,m), (x[(m+1):n]-x[1:(n-m)])/m))
  U = t(Htil) %*% t(Btil) %*% M %*% Btil %*% Htil
  V = U[-(1:m),-(1:m)]
  V = band(V,-(m-1),(m-1)) # Manually zero out elements
  return(Matrix(V,sparse=TRUE))
}

## Example functions

# Piecewise constant function
pieconst = function(x) {
  n = length(x)
  d = floor(n/5)
  f = rep(c(8,0,10,3,4), times=c(rep(d,4),n-4*d))
  return(f)
}

# Piecewise linear function
pielin = function(x) {
  x0 = seq(min(x),max(x),length=100)
  y0 = pielin100()
  return(approx(x0,y0,x,rule=2)$y)
}

pielin100 = function() {
  knots = matrix(0,6,2)
  knots[,1] = c(1,20,45,60,85,100)
  knots[,2] = c(1,2,6,8,5,6)
  f = rep(0,100)
  for (i in 1:(nrow(knots)-1)) {
    for (j in knots[i,1]:(knots[i+1,1])) {
      f[j] = (knots[i+1,1]-j)/(knots[i+1,1]-knots[i,1])*knots[i,2] +
        (j-knots[i,1])/(knots[i+1,1]-knots[i,1])*knots[i+1,2]
    }
  }

  f = 10*(f-min(f))/(max(f)-min(f))
  return(f)
}

# Piecewise quadratic function
piequad = function(x) {
  x0 = seq(min(x),max(x),length=100)
  y0 = piequad100()
  return(splinefun(x0,y0)(x))
}

piequad100 = function() {
  knots = matrix(0,4,2)
  knots[,1] = c(1,33,60,100)
  knots[,2] = c(8,6,2,4)
  f = rep(0,100); endval = 0
  for (i in 1:(nrow(knots)-1)) {
    sgn = (-1)^(i+1)
    mid = (knots[i,1]+knots[i+1,1])/2
    dif = knots[i+1,1]-knots[i,1]
    
    j = knots[i,1]
    intcp = endval - (sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
      (j-knots[i,1])/dif*knots[i+1,2])
    
    for (j in knots[i,1]:(knots[i+1,1])) {
      f[j] = intcp + sgn*n/5*(j-mid)^2/dif^2 + (knots[i+1,1]-j)/dif*knots[i,2] +
        (j-knots[i,1])/dif*knots[i+1,2]
    }
    endval = f[j]
  }
  
  f = 10*(f-min(f))/(max(f)-min(f))
  f = rev(max(f) - f + min(f))
  return(f)
}

# Smooth and wiggly function

smoothwiggly.fun = function(x) {
  f = function(a,b,c,x) return(a*(x-b)^2+c)
  fp = function(a,b,c,x) return(2*a*(x-b))

  a=-1; b=1/4; c=1;  
  if (x<=1/3) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=1/3;
  a=1; b=xx-fp(aa,bb,cc,xx)/(2*a); c=f(aa,bb,cc,xx)-a*(xx-b)^2;  
  if (x<=2/3) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=2/3;
  b=0.7; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.775) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.775;
  b=0.8; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.825) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.825;
  b=0.85; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.875) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.875;
  b=0.9; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.925) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.925;
  b=0.95; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.975) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.975;
  b=1; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  return(f(a,b,c,x))
}

smoothwiggly = function(x) {
  n = length(x); f = rep(0,n)
  for (i in 1:n) f[i] = smoothwiggly.fun(x[i])
  f = 10*(f-min(f))/(max(f)-min(f))
  return(f)
}
