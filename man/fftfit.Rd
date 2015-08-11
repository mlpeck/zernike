\name{fftfit}
\alias{fftfit}
\title{Fourier transform interferogram analysis}
\description{
  High level routines for FFT analysis of interferograms.
}
\usage{
fftfit(imagedata, cp = NULL, fringescale = 1, 
  sl = c(1, 1), filter = NULL, taper = 2, 
  zlist = makezlist(), zc0 = c(1:3, 6:7), 
  satarget = c(0, 0), astig.bath = c(0, 0), 
  puw.alg = "qual", uselm = FALSE, sgs = 3, plots = TRUE, CROP = FALSE)
}
\arguments{
  \item{imagedata}{A matrix containing the interferogram}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{fringescale}{Fringe spacing, in waves. Use 1 for single pass, 0.5 for double pass, etc.}
  \item{sl}{Position of sidelobe in the form c(x,y)}
  \item{filter}{Size of background filter around DC}
  \item{taper}{Size of taper applied to edge of half plane cut}
  \item{zlist}{Indexes of Zernike polynomials to fit to wavefront}
  \item{zc0}{Indexes of Zernike coefficients to be removed from net wavefront}
  \item{satarget}{Target 4th and 6th order SA coefficients in non-null tests of aspheres}
  \item{astig.bath}{Astigmatism coefficients for Bath geometry}
  \item{puw.alg}{Algorithm to use for phase unwrapping}
  \item{uselm}{Logical: use lm() for least squares fit}
  \item{sgs}{Sample Grid Spacing for least squares fits to wavefront values}
  \item{plots}{Logical: plot progress?}
  \item{CROP}{Center and crop maps?}
}
\details{
  If \code{is.null(filter)} (the default), \code{\link{pick.sidelobe}} will
  be called to select a Fourier domain sidelobe and background filter size.
  
  If \code{is.null(cp)} \code{\link{circle.pars}} is applied to the
  modulation to estimate the pupil parameters.
  
  See \code{\link{wf.net}} for details of the process of creating net and
  smoothed wavefronts from raw unwrapped wavefront maps.
  
  \code{puw.alg} Specifies the unwrapping algorithm. If \code{NULL} an algorithm
  that's likely to be successful will be selected. You can specify an algorithm
  by choosing \code{puw.alg=c("brcut", "ls", "lp", "modal", "qual")}.
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
\seealso{
  \code{\link{wf.net}},
  \code{\link{pupil.pars}},
  \code{\link{pick.sidelobe}}.
}
\keyword{mathematics}
