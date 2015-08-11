\name{readjpeg}
\alias{readjpeg}
\alias{readtiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a jpeg or tiff file}
\description{Reads a jpeg or tiff file and combines the channels
             to produce a monochrome image in a matrix.}
\usage{
readjpeg(filename, channels)
readtiff(filename, channels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{File name}
  \item{channels}{A vector of length 3 with the channel weights}
}
\details{Values in channels should be non-negative, but need not add to one.}
\value{A double precision matrix with the image data.}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{The matrix must have rows reversed and transposed to
   display properly with \code{image()}.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}
