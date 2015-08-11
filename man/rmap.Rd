\name{rmap}
\alias{rmap}
\alias{wrap}
\title{Utilities for phase unwrapping}
\description{
  Utility functions for use in 2D phase unwrapping.
}
\usage{
rmap(phase, dx = NULL, dy = NULL, plot = FALSE, ...)
wrap(phase)
}
\arguments{
  \item{phase}{Matrix of wrapped phase values}
  \item{dx}{Matrix of x differences}
  \item{dy}{Matrix of y differences}
  \item{plot}{Boolean: plot residue positions?}
  \item{...}{additional arguments for \code{image}}
}
\details{
  \code{dx} and \code{dy} must have the same dimension as \code{phase}.
}
\value{
  In \code{rmap} if \code{plot == TRUE}
  \item{nr}{the number of residues identified in the map}
  otherwise
  \item{phase}{wrapped phase returned by \code{wrap}}
  \item{residues}{Matrix the same size as \code{phase} with residues
    marked as + or - 1.}
}
\references{
  Ghiglia, D.C., and Pritt, M.D., 1998, \bold{Two-Dimensional Phase Unwrapping},
  New York: Wiley & Sons, Inc., ISBN 0-471-24935-1.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}.
  Steve Koehler is responsible for the efficient
  implementation of the \code{wrap} function.
}
\note{
	These are primarily intended for internal use
	but can be used interactively. Calling \code{rmap(phase, plot=TRUE)}
	will plot the positions of residues and return nothing.
	If \code{(plot==FALSE)} in the call to \code{rmap} a matrix
	the same size as \code{phase} is returned with residues identified
	with values of +1 or -1.
}
\seealso{
  Called by \code{\link{modalpuw}}, \code{\link{idiffpuw}}, \code{\link{brcutpuw}}.
}
\keyword{mathematics}
