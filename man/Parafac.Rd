% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Parafac.R
\name{Parafac}
\alias{Parafac}
\title{Performs parafac factorization via ALS}
\usage{
Parafac(Y, R)
}
\arguments{
\item{Y}{Array of dimension m1*m2*...*mK}

\item{R}{Desired rank of the factorization. Defalust is all columns(range=1:R)}
}
\value{
list with components
\item{U:}{List of basis vectors for parafac factorization, U[k]: mk*R for k=1,...,K. Columns of U[K] have norm 1 for k>1. }
\item{SqError:}{Vector of squared error for the approximation at each iteration (should be non-increasing)}
}
\description{
Performs parafac factorization via ALS
}
\examples{
A <- array(stats::rnorm(100*10*10), dim=c(100,10,10))
Parafac(A,4)
}
