\name{astig.bath}
\alias{astig.bath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zernike coefficients for astigmatism due to Bath astigmatism.
}
\description{
Calculates Bath astigmatism coefficients with optional rotation
of phi degrees.
}
\usage{
astig.bath(D, rc, s, lambda = 1e-06, phi = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{Diameter}
  \item{rc}{Radius of curvature}
  \item{s}{separation of reference and test beams}
  \item{lambda}{Wavelength -- defaults to 1 nm.}
  \item{phi}{angle of image horizontal relative to interferometer axis,
  in degrees}
}
\details{
  D, rc, s, and lambda must have the same units.
}
\value{
The Zernike coefficients for primary astigmatism terms.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\keyword{mathematics}
