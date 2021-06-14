\name{synth.interferogram}
\alias{synth.interferogram}
\title{Synthetic interferogram}
\description{
  Compute and display a synthetic interferogram.
}
\usage{
synth.interferogram(wf = NULL, zcoef = NULL, maxorder = NULL, 
  nr = nrow(wf), nc = ncol(wf), cp = NULL, 
  phi = 0, addzc = rep(0, 4), fringescale = 1, plots = TRUE)
}
\arguments{
  \item{wf}{A matrix of wavefront values}
  \item{zcoef}{A vector of Zernike coefficients}
  \item{maxorder}{Maximum Zernike polynomial order}
  \item{nr}{Number of rows in the output matrix}
  \item{nc}{Number of columns in the output matrix}
  \item{cp}{A list describing the pupil boundaries, as created by \code{\link{pupil.pars}}}
  \item{phi}{Amount to rotate the wavefront, in degrees}
  \item{addzc}{A 4-vector with piston, tilt, and defocus terms to be added}
  \item{fringescale}{Fringe scale. Should be 1 for single pass, 0.5 for double, etc.}
  \item{plots}{Logical: Plot the interferogram?}
}
\details{
  Either \code{wf} or \code{zcoef} should be non-null,
  but not both. If \code{zcoef} is specified \code{maxorder}
  must be as well.
  
  Additional piston, tilt, and defocus terms can be added
  to the calculated wavefront using \code{addzc}.
}
\value{
  A matrix of intensity levels in the calculated interferogram,
  assigned class \code{"\link{pupil}"}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{pupil}}.
}
\examples{
# create a list of zernikes
zcoef <- rnorm(length(zlist.fr$n), mean=0, sd=0.01)

iwf <- synth.interferogram(zcoef=zcoef, zlist=zlist.fr)

X11()

# show again with some tilt

iwf <- synth.interferogram(zcoef=zcoef, zlist=zlist.fr, addzc=c(0,5,5,0))
}
\keyword{graphics}
\keyword{hplot}
