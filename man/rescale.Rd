\name{rescale}
\alias{rescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rescale an image.}
\description{Rescale a matrix containing a bitmapped image
  using bilinear interpolation.}
\usage{
rescale(img, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{A matrix with image data.}
  \item{scale}{Scale factor.}
}
\details{A value <1 will shrink the image.}
\value{A matrix containing the rescaled image data.}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{\code{NA}'s are OK.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{Called by \code{\link{load.images}} if necessary.}

\keyword{array}
\keyword{arith}
