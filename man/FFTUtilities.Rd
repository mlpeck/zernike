\name{FFTUtilities}
\alias{wftophase}
\alias{padmatrix}
\alias{submatrix}
\alias{fftshift}
\title{FFT Utilities}
\description{
  Miscellaneous utilities for working with 
  2D images in the Fourier domain.
}
\usage{
wftophase(X, lambda=1)
padmatrix(X, npad, fill = mean(X, na.rm=TRUE))
submatrix(X, size = 255)
fftshift(X)
.up2(nr, nc=nr)
}
\arguments{
  \item{X}{A matrix}
  \item{lambda}{Value of the wavelength, in the same units as X}
  \item{npad}{Size of padded matrix}
  \item{fill}{Values to be assigned to padded matrix elements}
  \item{size}{Size of returned matrix}
  \item{nr}{A number}
  \item{nc}{A number}
}
\details{
  \code{wftophase} computes the complex phase from wavefront values.
  
  \code{padmatrix} pads a matrix to size \code{npad} x \code{npad},
  placing the original matrix in the lower left hand corner of the
  padded matrix.
  
  \code{submatrix} extracts a \code{size} x \code{size} matrix
  from the center of a larger matrix.
  
  \code{fftshift} shuffles the quadrants of a matrix around to
  put the DC element (1,1) in the center of the transformed matrix,
  with spatial frequencies increasing to the right and up.
}
\value{
  A matrix transformation of the input matrix \code{X}.
  
  \code{.up2} returns the next higher power of 2 than
  \code{max(nr, nc)}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These low level routines are used by several higher level
  functions that operate in the Fourier domain.
}
\seealso{
  \code{\link{startest}},
  \code{\link{fftfit}}.
}
\keyword{utilities}
