\name{pupil.rhotheta}
\alias{pupil.rhotheta}
\title{Polar coordinates}
\description{
  Calculate matrixes of polar coordinates for \code{\link{pupil}}'s.
}
\usage{
pupil.rhotheta(nrow, ncol, cp)
}
\arguments{
  \item{nrow}{Number of rows in interferogram images}
  \item{ncol}{Number of columns in interferogram images}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
}
\value{
  A list with the following components:
  \item{rho}{A matrix of radial coordinates}
  \item{theta}{A matrix of angular coordinates}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  My Zernike polynomial routines work in polar coordinates,
  which this function provides. Also, \code{NA}'s are used
  to fill the matrix outside the pupil boundary, making the
  returned values convenient for selecting pixels inside
  interferograms.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{pupil}}.
}
\keyword{utilities}
