\name{startest}
\alias{startest}
\title{Star test simulator}
\description{
  Simulates a star test.
}
\usage{
startest(wf=NULL, zcoef=NULL, zlist=makezlist(), phi=0,
	lambda = 1, defocus=5,
	nrow = 255, ncol = nrow, 
	cp = list(xc=128,yc=128,rx=127,ry=127,obstruct=0),
	obstruct=NULL, npad = 4, 
	gamma=2, psfmag=2, displaymtf=TRUE, displaywf=FALSE)
}
\arguments{
  \item{wf}{A matrix of class \code{pupil} containing wavefront values}
  \item{zcoef}{Vector of Zernike coefficients}
  \item{maxorder}{maximum Zernike polynomial order}
  \item{phi}{Angle to rotate wavefront}
  \item{lambda}{Wavelength, in same units as coefficients}
  \item{defocus}{Amount of defocus in waves}
  \item{nrow}{# rows in pupil matrix}
  \item{ncol}{# columns in pupil matrix}
  \item{cp}{pupil parameters}
  \item{obstruct}{Obstruction fraction}
  \item{npad}{Pad factor for FFT}
  \item{gamma}{Gamma value for graphics display}
  \item{psfmag}{Magnification factor for in focus PSF display}
  \item{displaymtf}{Logical: Display MTF?}
  \item{displaywf}{Logical: Display calculated wavefront?}
}
\details{
  If \code{wf} is \code{NULL} the wavefront is calculated from
  the the Zernike coefficients (which should be non-\code{NULL}).
}
\value{
  A list with the following components:
  \item{psf}{The in focus point spread function.}
  \item{otf}{The complex optical transfer function, a complex matrix of size \code{pupilsize}.}
  \item{mtf}{The modulation transfer function, a real matrix of size \code{pupilsize}.}
}
\references{  
  Born, M. and Wolf, E. 1999, \emph{Principles of Optics, 7th Edition},
  Cambridge University Press.

  Suiter, H. R., 1994, \emph{Star Testing Astronomical Telescopes},
  Willman-Bell, Inc.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{Zernike}},
  \code{\link{pupil}}.
}
\examples{
# a random, but probably almost diffraction limited, wavefront

temp <- startest(zcoef=rnorm(length(makezlist()$n), mean=0, sd=0.01), zlist=makezlist(), displaywf=TRUE)
}
\keyword{graphics}
\keyword{hplot}
