\name{zpm}
\alias{zpm}
\alias{zpm.arb}
\alias{zpmC}
\title{Matrixes of Zernike polynomials}
\description{
  Create a matrix of Zernike polynomial values.
}
\usage{
zpm(rho, theta, phi= 0 , maxorder = 14, nthreads=parallel::detectCores()/2)
zpmC(rho, theta, maxorder)
zapmC(rho, theta, maxorder)
zpmCP(rho, theta, maxorder)
zpm.arb(rho, theta, phi = 0, zlist = makezlist())
}
\arguments{
  \item{rho}{A vector of radial coordinates.}
  \item{theta}{A vector of angular coordinates, in radians.}
  \item{phi}{Orientation of the image, in degrees}
  \item{zlist}{A list of indexes, as returned by \code{\link{makezlist}}}
  \item{maxorder}{The maximum Zernike polynomial order}
  \item{nthreads}{Number of threads for threaded function call}
}
\details{
  \code{rho} and \code{theta} must be the same length.
}
\value{
  \code{zpm.arb} returns a matrix of size
  \code{length(rho)} x \code{length(zlist$n)} with values
  of Zernike polynomials
  evaluated at the polar coordinates
  \code{(rho, theta-pi*phi/180)}.
  
  \code{zpm}, \code{zpmC}, and \code{zpmCP} return a matrix of size
  \code{length(rho)} x \code{(maxorder/2+1)^2}
  of Zernike polynomial values including a piston term.
  
  \code{zapmC} returns a matrix of (approximate) Zernike annular polynomial values.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These are used by various routines to make least squares fits
  of sets of Zernike polynomials to measured wavefront values.

  \code{zpmC} is the C++ routine that does the computations
  for \code{zpm}. No column names are returned.
    
  Threaded computation of the matrix is now available using
  \code{zpmCP}.
  
}

\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{zlist.fr}},
  \code{\link{fitzernikes}}
}
\keyword{mathematics}

