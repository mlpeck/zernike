\name{psfit_options}
\alias{psfit_options}
\alias{get_zoptions}
\title{
Options for package \code{zernike}
}
\description{
Parameters controlling various aspects of PSI algorithms, Zernike polynomial fitting, and data display
}
\usage{
  psfit_options
  get_zoptions()
}
\format{
  Default parameter values are in the list \code{psfit_options}:
 \item{refine}     {TRUE}
 \item{puw_alg}    {"qual"}
 \item{fringescale}{1}
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
  \code{psfit_options} is defined in the package environment. To retrieve a copy use \code{get_zoptions}.
  The list can then be modified and passed to \code{\link{psifit}}.
  
  Parameters you might want to set include \code{zc0} which sets the target SA for ``numerical nulling.``
  
  \code{ptol} sets convergence tolerances for iterative PSI algorithms. These have different definitions and different
  values may be suitable for different algorithms.
  
  If you don't like the default color palette there are many other choices. If you like rainbows \code{\link{rygcb}}
  defined in this package produces a relatively perceptually uniform version.
  
  The number of variable Zernike terms in the algorith \code{\link{tiltpsi}} is controlled by \code{nzcs}. Set it to 3 to include defocus,
  5 to include primary astigmatism.
}
\author{M.L. Peck \email{mlpeck54@gmail.com}}
\keyword{list}
