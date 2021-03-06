\name{pupil.pars}
\alias{pupil.pars}
\title{Pupil parameters}
\description{
  Interactively determine the center, radius, and obstruction
  fraction of a circular or annular interferogram image.
}
\usage{
pupil.pars(im = NULL, obstructed = FALSE)
}
\arguments{
  \item{im}{A matrix containing an interferogram image}
  \item{obstructed}{Logical: is there a central obstruction?}
}
\details{
  In \code{pupil.pars}, if the image has already been plotted \code{im} can
  be \code{NULL}, which is the default.
}
\value{
  A list with the following components:
  \item{xc}{X coordinate of the center of the pupil}
  \item{yc}{Y coordinate of the center of the pupil}
  \item{rx}{Horizontal radius of the pupil}
  \item{ry}{Vertical radius of the pupil}
  \item{obstruct}{Obstruction fraction}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  \code{pupil.pars} uses the basic graphics library routine \code{\link{locator}}
  to interactively mark the edge of the pupil, and optionally
  the edge of the obstruction. After right clicking to terminate
  \code{locator()} a least squares fit is performed to the marked
  points to determine the center and radius of the pupil.
  
  Note that all routines that make use of Zernikes implicitly assume
  a circular pupil, or an annular one with small obstruction. We
  allow \code{rx != ry} for imaging sensors with non-square aspect
  ratios.
}
\seealso{
  Many routines require the pupil parameters in the form returned by
  \code{pupil.pars}. For example \code{\link{psifit}},
  \code{\link{fftfit}}, \code{\link{pupil}}, etc.
}
\keyword{utilities}
