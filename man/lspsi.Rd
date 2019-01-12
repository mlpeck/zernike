\name{lspsi}
\alias{lspsi}
\alias{lspsiC}
\title{Phase Shifting Interferometry}
\description{
  Least squares fitting of phase shifted interferograms.
}
\usage{
lspsi(images, phases, wt = rep(1, length(phases)))
lspsiC(images, phases, wt)
}
\arguments{
  \item{images}{An array containing the interferogram images}
  \item{phases}{A vector of phase shifts}
  \item{wt}{A vector of weights}
}
\details{
  \code{images} is a 3 dimensional array with dimensions \code{nrow x ncol x length(phases)},
  where \code{nrow} and \code{ncol} are the number of rows and columns in the individual
  interferogram images.
  
  \code{lspsi} reshapes the image array into a matrix and calls \code{lspsiC} which in turn calls the
  compiled C++ routine.
}
\value{
  A list containing the following components:
  \item{phi}{Estimated wrapped wavefront phase.}
  \item{mod}{Estimated modulation}
}
\references{
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
}
\seealso{
  \code{\link{psifit}}
}
\keyword{mathematics}
