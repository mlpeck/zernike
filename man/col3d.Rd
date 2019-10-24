\name{col3d}
\alias{col3d}
\title{OpenGL plot}
\description{
  Returns a vector of colors similar to \code{image()} display.
}
\usage{
col3d(surf, surf.col=topo.colors(256), zlim = NULL, eqa=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surf}{A matrix of surface values}
  \item{surf.col}{Color palette for surface}
  \item{zlim}{Range of values to display}
  \item{eqa}{Equal area per color?}
}
\value{
  A vector of color values the same length as \code{surf}.
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

\keyword{graphics}
\keyword{hplot}
