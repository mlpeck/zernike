\name{aiapsi}
\alias{aiapsi}
\alias{aiapsiC}
\alias{hkpsi}
\alias{tiltpsi}
\alias{tiltpsiC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Iterative algorithms for PSI with unknown phase shifts}
\description{Three iterative algorithms for PSI with unknown phase shifts.}
\usage{
aiapsi(im.mat, phases, ptol = 0.001, maxiter=20, trace=1)
aiapsiC(im.mat, phases_init, ptol, maxiter, trace)
hkpsi(im.mat, phases, maxiter = 20, ptol = 0.001,
  trace = 1, plotprogress = TRUE)
tiltpsi(im.mat, phases, coords, ptol = 0.01, maxiter = 20, trace = 1)
tiltpsiC(im.mat, phases, coords, ptol, maxiter, trace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im.mat}{a \emph{matrix} of interferogram values}
  \item{phases}{Starting guess for phase shifts}
  \item{ptol}{Convergence criterion for phase shifts}
  \item{maxiter}{Maximum number of iterations}
  \item{trace}{Boolean: Print some summary data at each iteration.}
  \item{plotprogress}{Plot some summary data for each iteration?}

  Also, for \code{tiltpsi} and \code{tiltpsiC}
  
  \item{coords}{Low order Zernike polynomial matrix}

}
\details{
  The ``variable tilt`` algorithm now allows an indefinite number of low order Zernike terms to be
  variable between phase steps. \code{coords} can be created with \code{\link{zpm}} setting \code{maxorder}
  to a small value, say 4, discarding the first (dc) column and retaining as many as desired. There must be at least
  two columns for tilts. The third will be defocus, the next two astigmatism, the next two primary coma, ...
  
  \code{aiapsi} and \code{tiltpsi} are wrappers for the calls to the C++ code in \code{aiapsiC} and \code{tiltpsiC}
  with sensible defaults for \code{ptol}, \code{maxiter}, and \code{trace}.
}
\value{
	A list containing the following elements:
	\item{phi}{The wrapped phase estimate. This is a vector as long as the number of rows in \code{im.mat}.}
	\item{mod}{Modulation estimate.}
	\item{phases}{Phase shift estimates.}
	\item{iter}{Number of iterations.}
	\item{sse}{Sum squared error at each iteration.}

	Also, for \code{tiltpsi}
	\item{zcs}{Matrix of Zernike coefficients, with one row for each column in \code{coords} and number of
	           columns = number of columns of im.mat.}
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
  \code{\link{psifit}}
}
\keyword{mathematics}
