% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dot_functions.R
\name{dot_functions}
\alias{dot_functions}
\alias{.divided_diff}
\alias{.b_mat_mult}
\alias{.h_mat_mult}
\title{In-place computations}
\usage{
.divided_diff(f, z)

.b_mat_mult(v, k, xd, tf_weighting, transpose, inverse)

.h_mat_mult(v, k, xd, di_weighting, transpose, inverse)
}
\arguments{
\item{f}{Function, or vector of function evaluations at the centers.}

\item{z}{Centers for the divided difference calculation.}

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

\item{di_weighting}{Should "discrete integration weighting" be used?
Multiplication by such a weighted H gives discrete integrals at the
design points; see details for more information. The default is \code{FALSE}.}
}
\value{
None. These functions \emph{overwrite} their input.
}
\description{
Each "dot" function accepts arguments as in its "non-dot"
counterpart, but peforms computations in place, overwriting the first input
argument (which must be a vector of the appropriate length) with the
desired output.
}
\details{
These functions should not be used unless you are intentionally
doing so for memory considerations and are nonetheless being very careful.

An \strong{important warning:} each "dot" function only works as expected if its
first argument is passed in as a vector of numeric type. If the first
argument is passed in as an integer vector, then since the output must (in
general) be a numeric vector, it cannot be computed in place (Rcpp performs
an implicit cast and copy when it converts this to NumericVector type for
use in C++).

Also, each "dot" function does not perform any error checking on its input
arguments. Use with care. More details on the computations performed by
individual functions are provided below.
}
\section{\code{.divided_diff()}}{

Overwrites \code{f} with all lower-order divided differences: each element \code{f[i]}
becomes the divided difference with respect to centers \code{z[1:i]}. See also
\code{\link[=divided_diff]{divided_diff()}}.
}

\section{\code{.b_mat_mult()}}{

Overwrites \code{v} with \code{B \%*\% v}, where \code{B} is the extended discrete derivative
matrix as returned by \code{b_mat()}. See also \code{\link[=b_mat_mult]{b_mat_mult()}}.
}

\section{\code{.h_mat_mult()}}{

Overwrites \code{v} with \code{H \%*\% v}, where \code{H} is the falling factorial basis
matrix as returned by \code{h_mat()}. See also \code{\link[=h_mat_mult]{h_mat_mult()}}.
}

\examples{
v = as.numeric(1:10) # Note: must be of numeric type
b_mat_mult(v, 1, 1:10) 
v
.b_mat_mult(v, 1, 1:10, FALSE, FALSE, FALSE)
v
}
