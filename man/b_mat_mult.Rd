% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_multiplication.R
\name{b_mat_mult}
\alias{b_mat_mult}
\title{Multiply by B matrix}
\usage{
b_mat_mult(v, k, xd, tf_weighting = FALSE, transpose = FALSE, inverse = FALSE)
}
\arguments{
\item{v}{Vector to be multiplied by B, the extended discrete derivative
matrix.}

\item{k}{Order for the extended discrete derivative matrix. Must be >= 0.}

\item{xd}{Design points. Must be sorted in increasing order, and have length
at least \code{k+1}.}

\item{tf_weighting}{Should "trend filtering weighting" be used? This is a
weighting of the discrete derivatives that is implicit in trend filtering;
see details for more information. The default is \code{FALSE}.}

\item{transpose}{Multiply by the transpose of B? The default is \code{FALSE}.}

\item{inverse}{Multiply by the inverse of B? The default is \code{FALSE}.}
}
\value{
Product of the extended discrete derivative matrix B and the input
vector \code{v}.
}
\description{
Multiplies a given vector by B, the extended discrete derivative matrix of a
given order, with respect to given design points.
}
\details{
The extended discrete derivative matrix of order \eqn{k}, with
respect to design points \eqn{x_1 < \ldots < x_n}, is denoted
\eqn{B^k_n}. It is square, having dimension \eqn{n \times n}. Acting on a
vector \eqn{v} of function evaluations at the design points, denoted \eqn{v
  = f(x_{1:n})}, it gives the discrete derivatives of \eqn{f} at the points
\eqn{x_{1:n}}:
\deqn{
  B^k_n v = (\Delta^k_n f) (x_{1:n}).
  }
The matrix \eqn{B^k_n} can be constructed recursively as the product of a
diagonally-weighted first difference matrix and \eqn{B^{k-1}_n}; see the
help file for \code{\link[=b_mat]{b_mat()}}, or Section 6.2 of Tibshirani (2020). Therefore,
multiplication by \eqn{B^k_n} or by its transpose can be performed in
\eqn{O(nk)} operations based on iterated weighted differences. See Appendix
D of Tibshirani (2020) for details.

The option \code{tf_weighting = TRUE} performs multiplication by \eqn{Z^k_n B^k_n}
where \eqn{Z^k_n} is an \eqn{n \times n} diagonal matrix whose top left
\eqn{k \times k} block equals the identity matrix and bottom right
\eqn{(n-k) \times (n-k)} block equals \eqn{W^k_n}, the latter being a
diagonal weight matrix that is implicit in trend filtering, as explained in
the help file for \code{\link[=d_mat_mult]{d_mat_mult()}}.

Lastly, the matrix \eqn{B^k_n} has a special \strong{inverse relationship} to the
falling factorial basis matrix \eqn{H^{k-1}_n} of degree \eqn{k-1} with
knots in \eqn{x_{k:(n-1)}}; it satisfies:
\deqn{
  Z^k_n B^k_n H^{k-1}_n = I_n,
  }
where \eqn{Z^k_n} is the \eqn{n \times n} diagonal matrix as described
above, and \eqn{I_n} is the \eqn{n \times n} identity matrix. This,
combined with the fact that the falling factorial basis matrix has an
efficient recursive representation in terms of weighted cumulative sums,
means that multiplying by \eqn{(B^k_n)^{-1}} or its transpose can be
performed in \eqn{O(nk)} operations. See Section 6.3 and Appendix D of
Tibshirani (2020) for details.
}
\examples{
v = sort(runif(10))
as.vector(b_mat(2, 1:10) \%*\% v)
b_mat_mult(v, 2, 1:10) 
}
\references{
Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.2.
}
\seealso{
\code{\link[=discrete_deriv]{discrete_deriv()}} for discrete differentiation at arbitrary query
points, \code{\link[=d_mat_mult]{d_mat_mult()}} for multiplying by the discrete derivative matrix,
and \code{\link[=b_mat]{b_mat()}} for constructing the extended discrete derivative matrix.
}
