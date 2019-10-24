\name{idiffpuw}
\alias{idiffpuw}
\title{Phase unwrapping by Integrating DIFFerences}
\description{
  Simple path following algorithm for two dimensional phase unwrapping.
}
\usage{
idiffpuw(phase, mask = phase, ucall = TRUE, dx = NULL, dy = NULL)
}
\arguments{
  \item{phase}{A matrix of wrapped phase values}
  \item{mask}{Matrix the same size as \code{phase} indicating masked pixels}
  \item{ucall}{Boolean: User call?}
  \item{dx}{Matrix of x differences}
  \item{dy}{Matrix of y differences}
}
\details{
  \code{mask} indicates pixels that shouldn't be unwrapped. In the simplest (default)
  case these are just pixels where phase is undefined.
}
\value{
  \code{if(ucall)}, a matrix of class ``\code{\link{pupil}}'' 
  with unwrapped wavefront values, otherwise a list with items:
    \item{puw}{Unwrapped phase}
    \item{uw}{Matrix indicating pixels that have been unwrapped.}
}
\references{
  Ghiglia, D.C., and Pritt, M.D., 1998, \bold{Two-Dimensional Phase Unwrapping},
  New York: Wiley & Sons, Inc., ISBN 0-471-24935-1.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}.
  Thanks to Steve Koehler for programming ideas to
  considerably speed up the algorithm.
}
\seealso{
  \code{\link{rmap}}, \code{\link{brcutpuw}}
}
\note{
  \code{\link{brcutpuw}} calls \code{\link{rmap}} first to
  check for the presence of residues. If there are none \code{idiffpuw}
  is guaranteed to work and is called to do the phase unwrapping.
  
  If there \emph{are} residues \code{brcutpuw} creates a mask then calls
  \code{idiffpuw} to unwrap unmasked portions of the phase map.
  
  This function is user callable as well; use a call
  of the form \code{idiffpuw(phase)}.
}
\keyword{mathematics}
