\name{circle.pars}
\alias{circle.pars}
\title{Pupil parameters}
\description{
  Automatically determine the center and radius of a circular interferogram image.
}
\usage{
circle.pars(im, fw=2, qt=0.995, excl=5, refine=2, 
    plots=TRUE, ask=TRUE, details=FALSE)
}
\arguments{
  \item{im}{A matrix containing an image of a circular disk}
  \item{fw}{Amount to smooth image}
  \item{qt}{Threshold to accept an edge point, expressed as a quantile}
  \item{excl}{number of pixels around border of frame to exclude}
  \item{refine}{radius range in pixels for a second pass estimate}
  \item{plots}{Plot edge candidates and fit?}
  \item{ask}{Wait for input before displaying fit?}
  \item{obstructed}{Logical: is there a central obstruction?}
}
\details{
  This routine partially implements the Canny algorithm for
  edge detection. After optionally smoothing the input image
  the gradient is calculated using a Sobel filter, and edge
  pixels are identified by locating local maxima in the
  magnitude of the gradient.
  
  The edge pixels with \code{qt} percentile largest gradients
  are passed to \code{lqs} in package \code{MASS} to determine
  robustly the best fit circle.
  
  Finally, if \code{refine > 0}, \emph{all} edge points
  within +- \code{refine} pixels of the previously determined
  edge are passed to \code{nls} for a second estimate of
  center point and radius.
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
  the smoothing parameter \code{fw} helps to suppress artifacts. Depending
  on how strong the actual edge is compared to artifacts \code{qt} may need
  to be either increased or decreased from the default value.
  
  if \code{details==TRUE} several more pieces of data are returned. This is
  mostly for debugging purposes and may be eliminated in the future.
}
\seealso{
  Many routines require the pupil parameters in the form returned by
  \code{circle.pars}. For example \code{\link{psifit}},
  \code{\link{fftfit}}, \code{\link{pupil}}, etc.
}
\keyword{utilities}
