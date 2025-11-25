\name{circle.auto}
\alias{circle.auto}
\title{Auto circular edge detection}
\description{
  Automatically determine the center and radius of a circular interferogram image.
}
\usage{
circle.auto(im, sigma=3, qt=0.995, excl=5,
    plots=TRUE, details=FALSE)
}
\arguments{
  \item{im}{A matrix containing an image of a circular disk}
  \item{sigma}{Amount to smooth image}
  \item{qt}{Threshold to accept an edge point, expressed as a quantile}
  \item{excl}{number of pixels around border of frame to exclude}
  \item{plots}{Plot edge candidates and fit?}
  \item{details}{Logical: return extra details?}
}
\details{
  This routine partially implements the Canny algorithm for
  edge detection. After optionally smoothing the input image
  the gradient is calculated using a Sobel filter, and edge
  pixels are identified by locating local maxima in the
  magnitude of the gradient.
  
  The edge pixels with \code{qt} percentile largest gradients
  are passed to \code{nlsrob} in package \code{robustbase} \code{lqs} in package \code{MASS} to determine
  robustly the best fit circle.
}
\value{
  A list with the following components:
  \item{xc}{X coordinate of the center of the pupil}
  \item{yc}{Y coordinate of the center of the pupil}
  \item{rx}{Horizontal radius of the pupil}
  \item{ry}{Vertical radius of the pupil = rx}
  \item{obstruct}{Obstruction fraction (always = 0)}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This routine is only effective on modulation estimates, and will almost
  certainly fail on interferogram images. Since data quality varies widely
  considerable experimentation may be needed on any given image. Increasing
  the smoothing parameter \code{sigma} helps to suppress artifacts. Depending
  on how strong the actual edge is compared to artifacts \code{qt} may need
  to be either increased or decreased from the default value.
  
  if \code{details==TRUE} several more pieces of data are returned. This is
  mostly for debugging purposes and may be eliminated in the future.
}
\seealso{
  Many routines require the pupil parameters in the form returned by
  \code{circle.auto}. For example \code{\link{psifit}},
  \code{\link{fftfit}}, \code{\link{pupil}}, etc.
}
\keyword{utilities}
