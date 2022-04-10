\name{crop}
\alias{crop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Crop an array}
\description{
  Crop a matrix or 3D array. Main application is to trim excess pixels
  from an image array, wavefront, etc.
}
\usage{
crop(img, cp, npad = 20, nxy=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{Array to be cropped.}
  \item{cp}{A list describing the pupil boundary.}
  \item{npad}{Amount of padding to leave around the edge.}
  \item{nxy}{Number pixels per side in cropped array.}
}
\details{
  \code{cp} is the list provided by \code{\link{circle.pars}}.
}
\value{
  \item{im}{The cropped array}
  \item{cp}{Revised value of cp}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\keyword{Utilities}
