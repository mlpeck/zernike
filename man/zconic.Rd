\name{zconic}
\alias{zconic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Zernike coefficients for a conic surface}
\description{
  Calculates the radially symmetric Zernike coefficient values up
  to order nmax for a conic surface relative to a sphere of the
  same paraxial radius of curvature.
}
\usage{
zconic(D, rc, b = -1, lambda = 1e-06, nmax = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{Diameter}
  \item{rc}{Radius of curvature}
  \item{b}{Conic constant}
  \item{lambda}{Wavelength -- defaults to 1 nm.}
  \item{nmax}{Maximum radial polynomial order}
}
\details{
  D, rc, and lambda must have the same units.
}
\value{
  A vector of length \code{nmax/2-1} of coefficient values, in increasing
  radial order, n=c(4,6, ...).
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{\code{\link{Zernike}}}
\examples{
	zconic(200,2000)
	zconic(10, 20, b=-1.05, lambda=632.8E-9, nmax=12)
}
\keyword{mathematics}% at least one, from doc/KEYWORDS
