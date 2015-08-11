\name{turbwf}
\alias{turbwf}
\title{Kolmogorov Turbulence}
\description{
  Simulates the optical effects of atmospheric turbulence using
  Noll's (1976) calculation of the covariance matrix of Zernike
  polynomials under Kolmogorov turbulence.
}
\usage{
turbwf(friedratio = 1, zlist = makezlist(2, 40), reps = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{friedratio}{Ratio of pupil diameter to Fried parameter}
  \item{zlist}{A list of Zernikes, as returned for example by \code{\link{makezlist}}}
  \item{reps}{Number of draws to simulate}
}
\details{
  The default value of zlist has 440 elements, which may be
  more than necessary for a reasonable representation of an
  ``atmospheric'' wavefront.
}
\value{
  A list with the following components:
  \item{zcoef.turb}{A \code{reps} x \code{length(zlist$n} matrix of simulated draws of Zernike coefficients.}
  \item{V}{Covariance matrix of the indexed Zernikes.}
 }
\references{
  Noll, R.J. 1976, \bold{Zernike polynomials and atmospheric turbulence}, \emph{J. Opt. Soc. Am.},
  Vol. 66, No. 3, p. 207.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  Typos in the original source material have been corrected in the code. Note that
  scintillation is not modelled.
}
\seealso{
  \code{\link{Zernike}},
  \code{\link{pupil}}.
}
\examples{
# Simulate a single draw from a turbulent atmosphere
zcoef.turb <- turbwf(friedratio=5, zlist=makezlist(2,30), reps=1)$zcoef.turb
# Warning: this can take a while
wf <- pupil(zcoef=zcoef.turb, zlist=makezlist(2,30))
plot(wf)
summary(wf)
}
\keyword{mathematics}
