\name{zpm}
\alias{zpm}
\alias{zpm.arb}
\alias{filldzm}
\alias{fillgradientzm}
\title{Matrixes of Zernike polynomials}
\description{
  Create a matrix of Zernike polynomial values, or their
  derivatives or gradient.
}
\usage{
zpm(rho, theta, phi= 0 , maxorder = 14)
zpm.arb(rho, theta, phi = 0, zlist = makezlist())
filldzm(rho, theta, phi = 0, zlist = makezlist())
fillgradientzm(rho, theta, phi = 0, zlist = makezlist())
}
\arguments{
  \item{rho}{A vector of radial coordinates.}
  \item{theta}{A vector of angular coordinates, in radians.}
  \item{phi}{Orientation of the image, in degrees}
  \item{zlist}{A list of indexes, as returned by \code{\link{makezlist}}}
  \item{maxorder}{The maximum Zernike polynomial order}
}
\details{
  \code{rho} and \code{theta} must be the same length.
}
\value{
  \code{zpm}, \code{zpm.arb} and \code{filldzm} return a matrix of size
  \code{length(rho)} x \code{length(zlist$n)} with values
  of Zernike polynomials or their radial derivatives
  evaluated at the polar coordinates
  \code{(rho, theta-pi*phi/180)}.
  
  \code{fillgradientzm} returns the gradient, in polar
  coordinates, of Zernikes in a \code{2*length(rho)} x
  \code{length(zlist$n)} matrix.
  Rows \code{1:length(rho)} contain the radial derivative,
  followed by \code{1/rho} times the tangential derivative.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These are used by various routines to make least squares fits
  of sets of Zernike polynomials to measured wavefront values.
}

\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.fr}},
  \code{\link{fitzernikes}}
}
\keyword{mathematics}

