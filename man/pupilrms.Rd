\name{pupilrms}
\alias{pupilrms}
\alias{pupilpv}
\alias{strehlratio}
\title{Wavefront statistics}
\description{
  Compute basic statistics of wavefronts stored in \code{"pupil"}
  objects.
}
\usage{
pupilrms(pupil)
pupilpv(pupil)
strehlratio(rms)
}
\arguments{
  \item{pupil}{A matrix of class \code{"pupil"}}
  \item{rms}{An rms wavefront error}
}
\value{
  Estimates of the RMS and P-V wavefront errors. \code{strehratio}
  calculates Mahajan's approximation to the Strehl ratio.
}
\references{
  Schroeder, D.J. 2000, \emph{Astronomical Optics, 2nd Edition},
  Academic Press, chapter 10.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  \code{pupilrms} simply calculates the standard deviation of
  finite values in the matrix \code{pupil}. This is a crude,
  but usually accurate enough estimate of the true RMS wavefront
  error.
  
  \code{\link{summary.pupil}} calls these functions.
}

\seealso{
  \code{\link{summary.pupil}}.
}
\examples{
zcoef <- rnorm(length(makezlist()$n), 0, 0.01)
wf <- pupil(zcoef=zcoef)
plot(wf)
summary(wf)
sqrt(crossprod(zcoef)) # A more accurate estimate of RMS
}
\keyword{statistics}
