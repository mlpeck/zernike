\name{fitzernikes}
\alias{fitzernikes}
\title{Least Squares fit to Zernike polynomials}
\description{
  Performs a least squares fit of a specified set of
  Zernike polynomials to a vector of wavefront measurements.
}
\usage{
fitzernikes(wf, rho, theta, eps=0, phi = 0, maxorder = 14, 
            nthreads = -1, isoseq = FALSE,
            usecirc = FALSE, ext_prec = FALSE)
}
\arguments{
  \item{wf}{A vector of wavefront values}
  \item{rho}{A vector of radial coordinates.}
  \item{theta}{A vector of angular coordinates, in radians.}
  \item{eps}{Obstruction fraction of annular aperture.}
  \item{phi}{Orientation of the image, in degrees}
  \item{maxorder}{Maximum Zernike polynomial order}
  \item{nthreads}{# threads for threaded function calls}
  \item{isoseq}{Boolean: use ISO/ANSI sequencing}
  \item{usecirc}{Boolean: Use circle Zernikes even if eps > 0?}
  \item{ext_prec}{Boolean: Use extended precision arithmetic for annular Zernike computation?}
}
\details{
  \code{wf}, \code{rho}, and \code{theta} must be the same length.
  
  As of version 3.7.0 Zernike polynomials in ISO/ANSI sequence can
  be used through a call to \code{\link{zpm_cart}}.
}
\value{
  The coefficients of the least squares fit.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  With the standard ordering of Zernikes Z0 is the piston
  term, Z1 and Z2 are x and y tilts, Z3 is defocus, etc.
}
\seealso{
  \code{\link{zpm}},
  \code{\link{zpm_cart}},
  \code{\link{zapm}},
  \code{\link{zapm_iso}},
  \code{\link{psifit}},
  \code{\link{fftfit}},
  \code{\link{vortexfit}},
  \code{\link{wf_net}}.
}
\keyword{statistics}
