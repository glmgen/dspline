% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_construction.R
\name{b_mat}
\alias{b_mat}
\title{Construct B matrix}
\usage{
b_mat(k, xd, tf_weighting = FALSE, row_idx = NULL)
}
\arguments{
\item{k}{Order for the extended discrete derivative matrix. Must be >= 0.}

\item{xd}{Design points. Must be sorted in increasing order, and have length
at least \code{k+1}.}

\item{tf_weighting}{Should "trend filtering weighting" be used? This is a
weighting of the discrete derivatives that is implicit in trend filtering;
see details for more information. The default is \code{FALSE}.}

\item{row_idx}{Vector of indices, a subset of \code{1:n} where \code{n = length(xd)},
that indicates which rows of the constructed matrix should be returned. The
default is \code{NULL}, which is taken to mean \code{1:n}.}
}
\value{
Sparse matrix of dimension \code{length(row_idx)} by \code{length(xd)}.
}
\description{
Constructs the extended discrete derivative matrix of a given order, with
respect to given design points.
}
\details{
The extended discrete derivative matrix of order \eqn{k}, with
respect to design points \eqn{x_1 < \ldots < x_n}, is denoted \eqn{B^k_n}.
It has dimension \eqn{n \times n}, and is banded with bandwidth \eqn{k+1}.
It can be constructed recursively, as follows. For \eqn{k \geq 1}, we first
define the \eqn{n \times n} extended difference matrix \eqn{\bar{B}_{n,k}}:
\deqn{
  \bar{B}_{n,k} =
  \left[\begin{array}{rrrrrrrrr}
  1 & 0 & \ldots & 0 & & & & \\
  0 & 1 & \ldots & 0 & & & & \\
  \vdots & & & & & & & \\
  0 & 0 & \ldots & 1 & & & & \\
  & & & -1 & 1 & 0 & \ldots & 0 & 0 \\
  & & & 0 & -1 & 1 & \ldots & 0 & 0 \\
  & & & \vdots & & & & & \\
  & & & 0 & 0 & 0 & \ldots & -1 & 1
  \end{array}\right]
  \begin{array}{ll}
  \left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
  \right\} & \hspace{-5pt} \text{$k$ rows} \\
  \left.\vphantom{\begin{array}{c} 1 \\ 0 \\ \vdots \\ 0 \end{array}}
  \right\} & \hspace{-5pt} \text{$n-k$ rows}
  \end{array}.
  }
We also define the \eqn{n \times n} extended diagonal weight matrix
\eqn{Z^k_n} to have first \eqn{k} diagonal entries equal to 1 and last
\eqn{n-k} diagonal entries equal to \eqn{(x_{i+k} - x_i) / k}, \eqn{i =
  1,\ldots,n-k}. The \eqn{k}th order extended discrete derivative matrix
\eqn{B^k_n} is then given by the recursion:
\deqn{
  \begin{aligned}
  B^1_n &= (Z^1_n)^{-1} \bar{B}_{n,1}, \\
  B^k_n &= (Z^k_n)^{-1} \bar{B}_{n,k} \, B^{k-1}_n,
  \quad \text{for $k \geq 2$}.
  \end{aligned}
  }
We note that the discrete derivative matrix \eqn{D^k_n} from \code{\link[=d_mat]{d_mat()}} is
simply given by the last \eqn{n-k} rows of the extended matrix \eqn{B^k_n}.

The option \code{tf_weighting = TRUE} returns \eqn{Z^k_n B^k_n} where \eqn{Z^k_n}
is the \eqn{n \times n} diagonal matrix as described above. This weighting
is implicit in trend filtering, as explained in the help file for
\code{\link[=d_mat_mult]{d_mat_mult()}}. See also Sections 6.1 and 6.2 of Tibshirani (2020) for
further discussion.

\strong{Note:} For multiplication of a given vector by \eqn{B^k_n}, instead of
forming \eqn{B^k_n} with the current function and then carrying out the
multiplication, one should instead use \code{\link[=b_mat_mult]{b_mat_mult()}}, as this will be more
efficient (both will be linear time, but the latter saves the cost of
forming any matrix in the first place).
}
\examples{
b_mat(2, 1:10)
b_mat(2, 1:10 / 10)
b_mat(2, 1:10, row_idx = 4:7)
}
\references{
Tibshirani (2020), "Divided differences, falling factorials, and
discrete splines: Another look at trend filtering and related problems",
Section 6.2.
}
\seealso{
\code{\link[=d_mat]{d_mat()}} for constructing the discrete derivative matrix, and
\code{\link[=b_mat_mult]{b_mat_mult()}} for multiplying by the extended discrete derivative matrix.
}
