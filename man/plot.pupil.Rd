\name{plot.pupil}
\alias{plot.pupil}
\alias{summary.pupil}
\title{Pupils and wavefronts}
\description{
  Plot and summary methods for objects of class \code{"pupil"}.
}
\usage{
plot.pupil(wf, cp=NULL, col = topo.colors(256), addContours = TRUE, cscale = FALSE,
   eqa=FALSE, zlim=NULL, ...)
summary.pupil(wf)
}
\arguments{
  \item{wf}{An object of class \code{"pupil"}}
  \item{cp}{Pupil parameters; a list as returned by \code{\link{pupil.pars}}}
  \item{col}{Color palette for plot}
  \item{addContours}{Logical: add contour lines?}
  \item{cscale}{Add a color scale legend?}
  \item{eqa}{Perform an "equal area" plot?}
  \item{zlim}{z limits to pass to image}
  \item{\dots}{Additional parameters to pass to \code{\link{image.default}}}
}
\details{
  These give simple plot and summary methods for objects of class \code{\link{pupil}}.
  
  If \code{eqa} is \code{TRUE}, each color in the palette will be used for an equal number
  of pixels (as opposed to representing an equal interval). Note: the color scale
  (when \code{cscale == TRUE}) may be inaccurate if a very small number of colors are used.
}
\value{
  none
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\seealso{
  \code{\link{pupil}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}},
  \code{\link{strehlratio}},
  \code{\link{pupil.pars}}.
}
\keyword{utilities}
