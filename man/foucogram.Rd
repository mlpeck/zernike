\name{foucogram}
\alias{foucogram}
\title{Simulate a Foucaultgram}
\description{
  Simulates the appearance of a wavefront under the
  Foucault test.
}
\usage{
foucogram(wf, edgex = 0, phradius = 0, slit = FALSE, 
  pad = 4, gamma = 1, map = FALSE, lev = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wf}{An object of class \code{\link{pupil}} containing wavefront values}
  \item{edgex}{lateral position of knife edge}
  \item{phradius}{radius of light source}
  \item{slit}{Logical: Is source a slit or pinhole?}
  \item{pad}{pad factor for FFT}
  \item{gamma}{Gamma value for graphics display}
  \item{map}{Logical: Overlay contours from wavefront map?}
  \item{lev}{Increment for contour levels, if used}
}
\details{
  The default value of 0 for \code{phradius} simulates a monochromatic point source.
  Try values in the range 10-30 to suppress diffraction effects.
}
\value{
  A matrix of intensity levels in the simulated image.
}
\references{
  See \url{http://home.netcom.com/~mpeck1/astro/foucault/ext_foucault.pdf} for
  an outline of the mathematical treatment of an extended source.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The key approximations here are treating the light source as monochromatic
  and spatially coherent, which is usually not the case for an extended source.
  Also, Fraunhofer diffraction theory is used.
}
\seealso{
  \code{\link{pupil}}
}
\keyword{graphics}
\keyword{hplot}
