\name{lspsi}
\alias{lspsi}
\title{Phase Shifting Interferometry}
\description{
  Least squares fitting of phase shifted interferograms.
}
\usage{
lspsi(images, phases, wt = rep(1, length(phases)))
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
}
\value{
  A list containing the following components:
  \item{phi}{Estimated wrapped wavefront phase; a matrix the same size as the interferogram images}
  \item{B}{Terms of the least squares fit; an \code{nrow x ncol x 3} array}
}
\references{
  Goldberg, K. and Bokor, J. 2001, \bold{Fourier-transform method
  of phase-shift determination}, \emph{Applied Optics},
  vol. 40, pp. 2886-2894; also available at
  \url{http://goldberg.lbl.gov/papers/Goldberg_AO_40(17).pdf}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  Appendix A of Goldberg and Bokor gives a useful summary of the
  least squares method of phase shifted interferometry analysis.
  The matrix B is defined in Equation A2 and solved in
  equation A5 of their paper.
}
\seealso{
  Called by \code{\link{psifit}},
  \code{\link{hkfit}}.
}
\keyword{mathematics}
