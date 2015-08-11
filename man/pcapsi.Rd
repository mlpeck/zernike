\name{pcapsi}
\alias{pcapsi}
\title{Vargas et al.'s Principal Components method for PSI
%%  ~~function to do ... ~~
}
\description{Compute the phase using the Principal components algorithm.}
\usage{
pcapsi(im.mat, BGSUB)
}
\arguments{
  \item{im.mat}{A \emph{matrix} of interferogram values}
  \item{BGSUB}{Logical - subtract the pixelwise mean as background estimate?}
}
\details{Images are input into an array by \code{\link{load.images}}. This must
be reshaped into a matrix for this function. Also, a mask should be applied if
available prior to the call.}
\value{
	A list containing the following elements:
	\item{phi}{The wrapped phase estimate. This is a vector as long as the number of rows in \code{im.mat}.}
	\item{mod}{Modulation estimate.}
	\item{phases}{Phase shift estimates.}
	\item{snr}{An estimate of the signal to noise ratio in the input data.}
	\item{eigen}{Singular values of the crossproduct matrix.}
}
\references{
J. Vargas, J. Antonio Quiroga, and T. Belenguer,
``Phase-shifting interferometry based on principal component analysis,''
\emph{Opt. Lett.} \bold{36}, 1326-1328 (2011)
\href{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-8-1326}{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-8-1326}

J. Vargas, J. Antonio Quiroga, and T. Belenguer,
``Analysis of the principal component algorithm in phase-shifting interferometry,''
\emph{Opt. Lett.} \bold{36}, 2215-2217 (2011)
\href{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-12-2215}{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-12-2215}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\seealso{
  \code{\link{pcafit}},
  \code{\link{wf.net}}
}
\keyword{mathematics}
