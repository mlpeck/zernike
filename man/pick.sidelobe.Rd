\name{pick.sidelobe}
\alias{pick.sidelobe}
\title{Select an interferogram sidelobe in the Fourier domain}
\description{
  Interactively locate the center of
  a first order sidelobe in the FFT of an interferogram,
  and mark the width of the background filter.
}
\usage{
pick.sidelobe(imagedata, logm=FALSE, gamma=3)
}
\arguments{
  \item{imagedata}{A matrix containing an interferogram image}
  \item{logm}{Logical: pass \code{fn="logMod"} to \code{\link{plot.cmat}}?}
  \item{gamma}{gamma value for display}
}
\details{
  Uses the basic graphics utility \code{\link{locator}}.
}
\value{
  A list with the following components:
  \item{sl}{The coordinates \code{c(x,y)} of the selected sidelobe}
  \item{filter}{Estimated size of background filter}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The high level FFT interferogram analysis routine \code{\link{fftfit}}
  requires the approximate location of the intended
  first order interferogram sidelobe to be specified.
}
\seealso{
  \code{\link{fftfit}},
}
\keyword{graphics}
\keyword{utilities}
