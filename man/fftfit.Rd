\name{fftfit}
\alias{fftfit}
\title{Fourier transform interferogram analysis}
\description{
  High level routines for FFT analysis of interferograms.
}
\usage{
fftfit(imagedata, cp = NULL, 
  sl = c(1, 1), filter = NULL, taper = 2, 
  options = psfit_options())
}
\arguments{
  \item{imagedata}{A matrix containing the interferogram}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{sl}{Position of sidelobe in the form c(x,y)}
  \item{filter}{Size of background filter around DC}
  \item{taper}{Size of taper applied to edge of half plane cut}
  \item{options}{a list of parameters passed to other functions. See \code{\link{psfit_options}}}.
}
\details{
  If \code{is.null(filter)} (the default), \code{\link{pick.sidelobe}} will
  be called to select a Fourier domain sidelobe and background filter size.
  
  If \code{is.null(cp)} \code{\link{circle.pars}} is applied to the
  modulation to estimate the pupil parameters.
  
  See \code{\link{wf_net}} for details of the process of creating net and
  smoothed wavefronts from raw unwrapped wavefront maps.
  
}
\value{
  A list with the following components:
  \item{phase}{Wrapped phase map}
  \item{mod}{The estimated modulation}
  \item{cp}{A list describing the pupil boundary}
  \item{cp.orig}{The precropped value of cp}
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
}
\references{
  Roddier, C. and Roddier, F. 1987, \bold{Interferogram
  analysis using Fourier transform techniques},
  \emph{Applied Optics}, vol. 26, pp. 1668-1673.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These functions are based largely on the work of Roddier and Roddier (1987).
}
\examples{
  # just run the example for 'vortexfit()'
  
  example(vortexfit, package="zernike", ask=FALSE)
}
\seealso{
  \code{\link{vortexfit}},
  \code{\link{wf_net}},
  \code{\link{pupil.pars}},
  \code{\link{pick.sidelobe}}.
}
\keyword{mathematics}
