\name{itfit}
\alias{itfit}
\title{Iterative algorithms for Phase Shifting Interferometry}
\description{
High level function that calls one of the iterative algorithms for PSI.
}
\usage{
itfit(images, phases, cp=NULL,
  maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE, REFINE=TRUE,
  tilts=NULL, tlim=0.5, nlpref=1,
  it.alg = c("aia", "hk", "tilt"),
  fringescale=1, zlist=makezlist(), zc0=c(1:3, 6:7), 
  satarget=c(0,0), astig.bath=c(0,0)
  puw.alg = "qual", uselm=FALSE, sgs=1, plots=TRUE, CROP=FALSE)
}
\arguments{
  \item{images}{An array containing the interferogram images}
  \item{phases}{A vector of phase shifts}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{maxiter}{maximum number of iterations of phase shift updates}
  \item{ptol}{A measure of phase shift convergence tolerance}
  \item{trace}{print trace of phase step estimates}
  \item{plotprogress}{Logical - plot progress of phase shift estimation?}
  \item{REFINE}{Logical - run the algorithm a second time after calculating mask}
  \item{tilts}{initial guess of tilts}
  \item{tlim}{tilt value limit}
  \item{nlpref}{preferred nonlinear optimizer}
  \item{it.alg}{iterative algorithm to use}
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
  interferogram mask automatically. If \code{REFINE == TRUE} the selected iterative
  algorithm will be called a second time with the masked data.
  Selecting \code{cp = NULL} will generate an error if \code{it.alg == "tilt"}.

  The arguments \code{tilts, tlim} and \code{nlpref} are passed to \code{\link{tiltpsi}}
  if \code{it.alg == "tilt"} and are ignored otherwise.

  \code{puw.alg} specifies the unwrapping algorithm. You can specify an algorithm
  by choosing \code{puw.alg=c("brcut", "ls", "lp", "modal", "qual")}.
  
  See \code{\link{wf.net}} for details of the process of creating net and
  smoothed wavefronts from raw unwrapped wavefront maps.
  
}
\value{
  A list with the following components:
  \item{phase}{Raw, wrapped phase map}
  \item{mod}{The estimated modulation}
  \item{phases}{final phase shift estimates}
  \item{cp}{A list describing the pupil boundary}
  \item{cp.orig}{Uncropped value of cp}
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
  \item{iter}{Number of iterations of the algorithm}
  \item{sse}{vector of sum squared errors}

  And, if \code{it.alg == "tilt"}
  \item{tilts}{final estimate of tilts for each frame}
}
\references{
Zhaoyang Wang and Bongtae Han,
``Advanced iterative algorithm for phase extraction of randomly phase-shifted interferograms,''
\emph{Opt. Lett.} 29, 1671-1673 (2004).

Han, G-S and Kim, S-W,, ``Numerical correction of reference phases
in phase-shifting interferometry by iterative least squares fitting,''
\emph{Applied Optics} 33, 7321-7325 (1994),

Lin, B-J et al., ``An iterative tilt-immune phase-shifting algorithm,''
OSA conference Optical Fabrication and Testing 2010.

}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\seealso{
  \code{\link{psifit}},
  \code{\link{pcafit}},
  \code{\link{wf.net}},
  \code{\link{circle.pars}}.
}
\keyword{mathematics}
