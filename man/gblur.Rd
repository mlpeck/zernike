\name{gblur}
\alias{gblur}
\title{Gaussian blur}
\description{
  Blur an image by sigma pixels
}
\usage{
gblur(X, sigma)
gblur_complex(X, sigma)
}
\arguments{
  \item{X}{A matrix representing an image, or a complex valued matrix}
  \item{sigma}{Width of the Gaussian convolution kernel, in pixels}
}
\value{
  The filtered matrix X.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{convolve2d}}
}
\keyword{image processing}
