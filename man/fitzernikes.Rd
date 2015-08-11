\name{fitzernikes}
\alias{fitzernikes}
\title{Least Squares fit to Zernike polynomials}
\description{
  Performs a least squares fit of a specified set of
  Zernike polynomials to a vector of wavefront measurements.
}
\usage{
fitzernikes(wf, rho, theta, phi = 0, maxorder = 14, uselm = FALSE)
}
\arguments{
  \item{wf}{A vector of wavefront values}
  \item{rho}{A vector of radial coordinates.}
  \item{theta}{A vector of angular coordinates, in radians.}
  \item{phi}{Orientation of the image, in degrees}
  \item{maxorder}{Maximum Zernike polynomial order}
  \item{uselm}{Boolean: use lm() for least squares fit}
}
\details{
  \code{wf}, \code{rho}, and \code{theta} must be the same length.
}
\value{
  The model fit as returned by \code{\link{lm}},
  or the coefficients of the least squares fit if
  \code{uselm} is \code{FALSE}.  
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The model fit is of the form \code{wf ~ Z0+Z1+Z2+\dots}.
  With the standard ordering of Zernikes Z0 is the piston
  term, Z1 and Z2 are x and y tilts, Z3 is defocus, etc.
}
\seealso{
  \code{\link{zpm}},
  \code{\link{psifit}},
  \code{\link{fftfit}},
  \code{\link{wf.net}}.
}
\keyword{statistics}
