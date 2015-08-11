\name{modalfit}
\alias{modalfit}
\alias{zonalfit}
\title{Foucault test data reduction}
\description{
  Data reduction for moving source Foucault tests.
}
\usage{
modalfit(rho, f, rc, R, run = NULL, theta = rep(0, length(rho)), 
  zlist = list(n = seq(4, 12, by = 2), m = rep(0, 5), t = rep("n", 5)))
zonalfit(rho, f, rc, R, 
  zlist = list(n = seq(4, 12, by = 2), m = rep(0, 5), t = rep("n", 5)))
}
\arguments{
  \item{rho}{A vector of relative zonal radii}
  \item{f}{A vector of longitudinal aberration measurements}
  \item{rc}{Radius of curvature of the mirror}
  \item{R}{Semi-diameter of the mirror}
  \item{run}{A vector of factor levels, if different measurement runs are included}
  \item{theta}{Orientation angle of the mirror, if different orientations are present}
  \item{zlist}{Indexes of Zernikes to be fit}
}
\details{
  All linear measurements must be in millimeters.
}
\value{
  A long list describing the fit to the data. Please contact the
  author if you need details.
}
\references{
  Modal analysis of Foucault test data is described by the author in
  \url{http://home.netcom.com/~mpeck1/astro/modal/modal.pdf}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These are provided for completeness only. Please contact me if you
  need details on their use.
}
\keyword{statistics}

