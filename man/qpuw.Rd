\name{qpuw}
\alias{qpuw}
\title{Quality guided algorithm for phase unwrapping}
\description{
  Quality guided algorithm for two dimensional phase unwrapping.
}
\usage{
  qpuw(phase, qual)
}
\arguments{
  \item{phase}{A matrix of wrapped phase values}
  \item{qual}{A matrix of quality values the same size as \code{phase}.}
}
\value{
  \item{puw}{A matrix of class \code{"\link{pupil}"} with the unwrapped wavefront.} 
}
\references{
  Ghiglia, D.C., and Pritt, M.D., 1998, \bold{Two-Dimensional Phase Unwrapping},
  New York: Wiley & Sons, Inc., ISBN 0-471-24935-1.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\seealso{\code{\link{idiffpuw}}, \code{\link{brcutpuw}}} 
\note{
  This is a straightforward implementation of the quality guided algorithm
  of G&P.
}
\keyword{mathematics}

