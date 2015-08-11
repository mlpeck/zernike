\name{convolve2d}
\alias{convolve2d}
\title{2D convolution}
\description{
  General 2D convolution using FFTs
}
\usage{
convolve2d(im, kern)
}
\arguments{
  \item{im}{A matrix representing an image}
  \item{kern}{the convolution kernel}
}
\details{
}
\value{
  The filtered matrix im.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{gblur}}. Called by \code{\link{circle.pars}}.
}
\keyword{mathematics}
