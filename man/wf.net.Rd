\name{wf.net}
\alias{wf.net}
\title{Wavefront smoothing}
\description{
  Calculate net and smoothed wavefronts from a raw wavefront containing
  low order nuisance aberrations.
}
\usage{
wf.net(wf.raw, rho=NULL, theta=NULL, cp, sgs=3, zlist=makezlist(),
    zc0=c(1:3, 6:7), satarget=c(0,0), astig.bath=c(0,0),
    uselm=FALSE, plots=TRUE)
}
\arguments{
  \item{wf.raw}{Raw wavefront to be processed}
  \item{rho}{Radial coordinates in pupil as returned by \code{\link{pupil.rhotheta}}}
  \item{theta}{Angular coordinates in pupil as returned by \code{\link{pupil.rhotheta}}}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{sgs}{Sample Grid Spacing for least squares fits to wavefront values}
  \item{zlist}{Indexes of Zernike polynomials to fit to wavefront}
  \item{zc0}{Indexes of Zernike coefficients to be removed from net wavefront}
  \item{satarget}{Target 4th and 6th order SA coefficients in non-null tests of aspheres}
  \item{astig.bath}{Astigmatism coefficients for Bath geometry}
  \item{uselm}{Boolean: use lm() for least squares fit}
  \item{plots}{Logical: plot calculated wavefronts?}
}
\details{
  In performing least squares fits to the raw data the pupil
  is sampled on a square grid with spacing
  \code{sgs} by \code{sgs} pixels. The default value of 3 appears to
  give a sufficiently dense sample for typical video resolution images
  without using too much memory or CPU cycles.
  
  Passing \code{rho} and \code{theta} is optional. If \code{rho}
  is \code{NULL} they are calculated using \code{\link{pupil.rhotheta}}.
  
  In a non-null test of an asphere setting \code{satarget} to the desired Zernike
  coefficients for 3rd and 5th order spherical aberration will produce a net wavefront
  with the target SA coefficients removed.
}
\value{
  A list with the following components:
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\seealso{
  Called by
  \code{\link{fftfit}},
  \code{\link{psifit}},
  \code{\link{hkfit}},
  \code{\link{pcafit}},
  \code{\link{aiafit}}.
  
  Calls
  \code{\link{fitzernikes}},
  \code{\link{pupil}}.
}
\keyword{mathematics}
\keyword{graphics}
