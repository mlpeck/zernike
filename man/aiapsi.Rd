\name{aiapsi}
\alias{aiapsi}
\alias{hkpsi}
\alias{tiltpsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Iterative algorithms for PSI with unknown phase shifts}
\description{Three iterative algorithms for PSI with unknown phase shifts.}
\usage{
aiapsi(im.mat, phases, maxiter = 20, ptol = 0.001,
  trace = 1, plotprogress = TRUE)
hkpsi(im.mat, phases, maxiter = 20, ptol = 0.001,
  trace = 1, plotprogress = TRUE)
tiltpsi(im.mat, phases, x, y, tilts = NULL, nlpref = 1, tlim = 0.5,
  maxiter = 20, ptol = 0.001, trace = 1, plotprogress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im.mat}{a \emph{matrix} of interferogram values}
  \item{phases}{Starting guess for phase shifts}
  \item{maxiter}{Maximum number of iterations}
  \item{ptol}{Convergence criterion for phase shifts}
  \item{trace}{Print some summary data every trace'th iteration.}
  \item{plotprogress}{Plot some summary data for each iteration?}

  Also, for \code{tiltpsi}

  \item{x}{x coordinate for each pixel.}
  \item{y}{y coordinate for each pixel.}
  \item{tilts}{Starting guess for tilts.}
  \item{nlpref}{Preferred nonlinear optimizer.}
  \item{tlim}{Maximum tilt difference from overall tilt values.}
}
\details{
The default of \code{tilts = NULL} will set the starting guess
for all tilts to 0.

\code{nlpref} picks the preferred nonlinear optimizer in function \code{tiltpsi}.
The default choice uses the base package function \code{nlminb}. If \code{nlpref > 1}
the function \code{solnp} from package \code{Rsolnp} is used instead.
}
\value{
	A list containing the following elements:
	\item{phi}{The wrapped phase estimate. This is a vector as long as the number of rows in \code{im.mat}.}
	\item{mod}{Modulation estimate.}
	\item{phases}{Phase shift estimates.}
	\item{iter}{Number of iterations.}
	\item{sse}{Sum squared error at each iteration.}

	Also, for \code{tiltpsi}
	\item{tilts}{Final tilt estimates, given as offsets from frame 1 tilts.}
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
  \code{\link{itfit}}
}
\keyword{mathematics}
