\name{psfit_options}
\alias{psfit_options}
\title{
Options for package \code{zernike}
}
\description{
Get and optionally set parameters controlling various aspects of PSI algorithms, Zernike polynomial fitting, and data display
}
\usage{
  psfit_options(...)
}
\arguments{
 \item{refine}{TRUE}
 \item{puw_alg}{"qual"}
 \item{fringescale}{1}
 \item{wt}{NULL}
 \item{bgsub}{TRUE}
 \item{maxiter    }{20}
 \item{ptol       }{1e-04}
 \item{trace      }{1}
 \item{nzcs       }{2}
 \item{zc0        }{c(1, 2, 3, 6, 7)}
 \item{satarget   }{c(0, 0)}
 \item{astig.bath }{c(0, 0)}
 \item{maxorder   }{14}
 \item{uselm      }{FALSE}
 \item{sgs        }{1}
 \item{plots      }{TRUE}
 \item{crop       }{FALSE}
 \item{colors     }{topo.colors(256)}
}
\details{
  Calling \code{psfit_options} with an empty argument list returns the default values of the options used in \code{\link{psifit}}.
  The list can be modified directly or by passing argument value pairs to the function call.
  
  Parameters you might want to change include:
  
  \code{zc0} sets the target SA for ``numerical nulling.`` This is a vector of length 2 setting the target values of primary and 5th order SA.
  
  \code{ptol} sets convergence tolerances for iterative PSI algorithms. These have different definitions and different
  values may be suitable for different algorithms. A value around 0.01 is appropriate for \code{\link{tiltpsi}}.
  
  The number of variable Zernike terms in the algorithm \code{\link{tiltpsi}} is controlled by \code{nzcs}. Set it to 3 to include defocus,
  5 to include primary astigmatism, 7 to include coma.
  
  \code{maxorder} sets the maximum Zernike polynomial order for wavefront fitting. It must be even and at least 6. The default generally produces a
  good wavefront representation but you may want to experiment with higher order fits.

  If you don't like the default color palette there are many other choices. If you like rainbows \code{\link{rygcb}}
  defined in this package produces a relatively perceptually uniform version that's well suited for display on an RGB monitor.
  
}
\value{
  A named list with the current values of the arguments.
}
\author{M.L. Peck \email{mlpeck54@gmail.com}}
\keyword{list}
