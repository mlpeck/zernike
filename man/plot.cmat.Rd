\name{plot.cmat}
\alias{plot.cmat}
\title{Plot a complex matrix}
\description{
  Plot a real valued function of a complex matrix
}
\usage{
plot.cmat(X, fn = "Mod", col = grey256, 
	cp=NULL, zoom=1, gamma=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A complex valued matrix}
  \item{fn}{A function returning a real value}
  \item{col}{Color palette for graph}
  \item{cp}{pupil parameters as returned by \code{\link{pupil.pars}}}
  \item{zoom}{zoom factor for display}
  \item{gamma}{gamma value for display}
  \item{\dots}{Other parameters to pass to \code{\link{image.default}}}
}
\details{
  In addition to the functions described in \code{\link{complex}}
  \code{fn} can be assigned the values \code{"logMod"}, which
  will call an internally defined function returning the value
  \code{log(1+Mod(X))}, \code{"Mod2"} to plot the power spectrum,
  and \code{"logMod2"} to plot the logarithm of the power spectrum.
  
  If the parameter cp is passed axes will display spatial frequencies
  in cycles per pupil radius.
}
\value{
  none
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  This is used primarily for displaying FFT's of
  interferograms. In the case of an interferogram
  in which the background has not been removed use
  \code{fn="logMod"} to make the first order sidelobes
  visible.
}
\seealso{
  \code{\link{pick.sidelobe}},
  \code{\link{fftfit}}.
}
\keyword{graphics}
