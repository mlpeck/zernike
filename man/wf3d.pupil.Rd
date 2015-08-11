\name{wf3d.pupil}
\alias{wf3d}
\title{OpenGL wavefront plot}
\description{
  Interactive plot of a wavefront using the OpenGL package
  \pkg{rgl}. This is a 3D plotting method for objects of
  class \code{"\link{pupil}"}.
}
\usage{
wf3d.pupil(wf, cp=NULL, zoom.wf = 1, surf.col = topo.colors(256), bg.col = "black",
    eqa=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wf}{A matrix of wavefront values}
  \item{cp}{A list describing the pupil boundary}
  \item{zoom.wf}{Zoom factor for heights}
  \item{surf.col}{Color palette for surface}
  \item{bg.col}{Background color}
  \item{eqa}{Equal area per color?}
}
\details{
  The default color palette will match the colors in the default
  version of \code{\link{plot.pupil}}.
}
\value{
  none
}
\references{
  The \pkg{rgl} package is described at 
  \url{http://rgl.neoscientists.org/about.shtml},
  and available from CRAN.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{
  \code{\link{plot.pupil}}
}
\examples{
# create a random wavefront

wf <- pupil(zcoef=rnorm(length(makezlist()$n), mean=0, sd=0.01))
# the default method

plot(wf)

#this is more fun

wf3d(wf)
}

\keyword{graphics}
\keyword{hplot}
