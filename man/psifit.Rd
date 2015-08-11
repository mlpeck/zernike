\name{psifit}
\alias{psifit}
\title{Phase Shifting Interferometry}
\description{
  High level function for Least squares analysis of
  phase shifted interferograms.
}
\usage{
psifit(images, phases, cp=NULL, wt=rep(1,length(phases)), 
  fringescale=1, zlist=makezlist(), 
  zc0=c(1:3, 6:7), satarget=c(0,0), astig.bath=c(0,0)
  puw.alg = "qual", uselm=FALSE, sgs=2, plots=TRUE, CROP=FALSE)
}
\arguments{
  \item{images}{An array containing the interferogram images}
  \item{phases}{A vector of phase shifts}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{wt}{\code{psifit} only: A vector of weights}
  \item{fringescale}{Fringe spacing, in waves. Use 1 for single pass, 0.5 for double pass, etc.}
  \item{zlist}{Indexes of Zernike polynomials to fit to wavefront}
  \item{zc0}{Indexes of Zernike coefficients to be removed from net wavefront}
  \item{satarget}{Target 4th and 6th order SA coefficients in non-null tests of aspheres}
  \item{astig.bath}{Astigmatism coefficients for Bath geometry}
  \item{puw.alg}{Algorithm to use for phase unwrapping}
  \item{uselm}{Boolean: use lm() for least squares fit}
  \item{sgs}{Sample Grid Spacing for least squares fits to wavefront values}
  \item{plots}{Logical: plot progress?}
  \item{CROP}{Center and crop maps?}
}
\details{
  \code{images} is a 3 dimensional array with dimensions \code{nrow x ncol x length(phases)},
  where \code{nrow} and \code{ncol} are the number of rows and columns in the individual
  interferogram images.
  
  If \code{cp == NULL} \code{\link{circle.pars}} is called to construct the
  interferogram mask automatically.
  
  See \code{\link{wf.net}} for details of the process of creating net and
  smoothed wavefronts from raw unwrapped wavefront maps.
  
  \code{puw.alg} specifies the unwrapping algorithm. You can specify an algorithm
  by choosing
  \code{puw.alg=c("brcut", "ls", "lp", "modal", "qual")}.

}
\value{
  A list with the following components:
  \item{phase}{Raw, wrapped phase map}
  \item{mod}{The estimated modulation}
  \item{cp}{A list describing the pupil boundary}
  \item{cp.orig}{Uncropped value of cp}
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
The low level function implementing the least squares algorithm
is \code{\link{lspsi}}.
}

\seealso{
  \code{\link{lspsi}},
  \code{\link{wf.net}},
  \code{\link{hkfit}},
  \code{\link{aiafit}},
  \code{\link{pcafit}}.
}
\keyword{mathematics}
