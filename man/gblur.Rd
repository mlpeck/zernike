\name{gblur}
\alias{gblur}
\title{Gaussian blur}
\description{
  Blur an image by fw pixels
}
\usage{
gblur(X, fw=0, details=FALSE)
}
\arguments{
  \item{X}{A matrix representing an image}
  \item{fw}{Width of the Gaussian convolution kernel, in pixels}
  \item{details}{Return convolution kernel?}
}
\details{
  \code{fw} is the standard deviation of the Gaussian.
}
\value{
  The filtered matrix X.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  the \code{details} option is mostly for debugging purposes and may go away.
}
\seealso{
  \code{\link{convolve2d}}
}
\keyword{mathematics}
