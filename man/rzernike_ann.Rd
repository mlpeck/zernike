% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rzernike_ann}
\alias{rzernike_ann}
\title{Radial Zernike Annular polynomials}
\usage{
rzernike_ann(rho, eps, n, m, xq, qwts)
}
\arguments{
\item{rho}{a vector of radial coordinates.}

\item{eps}{the obstruction fraction 0 <= eps < 1.}

\item{n}{the maximum radial order required}

\item{m}{azimuthal order}

\item{xq}{nodes for Gauss-Legendre quadrature}

\item{qwts}{weights for Gauss-Legendre quadrature}
}
\value{
A length(rho) x (n-m)/2+1 column matrix of radial Zernike Annular polynomial values evaluated at the input
radial coordinates. The radial indexes are in increasing order from m, m+2, ..., n.
}
\description{
Create a matrix of Zernike Annular polynomial values in
extended Fringe sequence for a set of polar coordinates.
}
\details{
To the author's knowledge no recurrence relations for radial Zernike annular polynomials
have been published, even though several are well known for the closely related Zernike circle polynomials.
However the m=0 polynomials representing axially symmetric aberrations are just shifted Legendre polynomials
with an easily derived recurrence relation. This routine makes use of that fact to generate
recurrence relations for arbitrary polynomial indexes using chebyshev's algorithm with modified moments.
The modified moments are calculated using Gauss-Legendre quadrature. If enough quadrature nodes
were chosen the quadrature is nominally exact, as are the resulting annular Zernike values.
}
\seealso{
This function is called by \code{\link[=zapm]{zapm()}} and \code{\link[=zapm_iso]{zapm_iso()}}.
}
