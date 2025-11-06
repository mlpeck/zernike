\name{zernike_options}
\alias{zernike_options}
\title{
Options for routines in package `zernike'.
}
\description{
Get and optionally set parameters controlling various aspects of PSI algorithms, Zernike polynomial fitting, and data display
}
\usage{
  zernike_options(...)
}
\arguments{
 \item{colors     }{topo.colors(256)}
 \item{plots      }{TRUE}
 \item{refine     }{TRUE}
 \item{puw_alg    }{"qual"}
 \item{fringescale}{1}
 \item{wt         }{NULL}
 \item{bgsub      }{TRUE}
 \item{maxiter    }{20}
 \item{ptol       }{1e-04}
 \item{trace      }{1}
 \item{nzcs       }{2}
 \item{zc0        }{6:7}
 \item{satarget   }{c(0, 0)}
 \item{astig.bath }{c(0, 0)}
 \item{maxorder   }{14}
 \item{isoseq     }{FALSE}
 \item{usecirc    }{FALSE}
 \item{nthreads   }{parallel::detectCores()/2}
 \item{crop       }{FALSE}
 \item{sigma_ed   }{3.}
 \item{qt_ed      }{0.995}
}
\details{
  Calling \code{zernike_options} with an empty argument list returns the default values of the options used in \code{\link{psifit}}
  and \code{\link{wf_net}} as itemized above.
  The list can be modified directly or by passing argument value pairs to the function call.
  
  Parameters you might want to change include:
  
  \code{satarget} sets the target SA for ``numerical nulling.`` This is a vector of length 2 setting the target values of primary and 5th order SA.
  
  \code{ptol} sets convergence tolerances for iterative PSI algorithms. These have different definitions and different
  values may be suitable for different algorithms. A value around 0.01 is appropriate for \code{\link{tiltpsi}}.
  
  The number of variable Zernike terms in the algorithm \code{\link{tiltpsi}} is controlled by \code{nzcs}. Set it to 3 to include defocus,
  5 to include primary astigmatism, 7 to include coma.
  
  \code{maxorder} sets the maximum Zernike polynomial order for wavefront fitting. It must be even and at least 6. The default generally produces a
  good wavefront representation but you may want to experiment with higher order fits.

  Parameters \code{isoseq} and \code{usecirc} control which sets of Zernike (annular) polynomials are called for wavefront fitting. If Zernikes in ANSI/ISO sequence
  are desired set \code{isoseq = TRUE}. The default \code{FALSE} uses extended Fringe ordering. Set \code{usecirc = TRUE} if you want to use Zernike circle polynomials
  even for obstructed wavefronts.
  
  \code{puw_alg} specifies the phase unwrapping algorithm. Current options are "qual", "brcut", "lpbrcut", and "lp".
  
  A new and somewhat experimental feature from version 3.6.0 of the package is threaded computation of Zernike polynomial matrices. The number of threads
  used is set with the option \code{nthreads}. The default is to use half the number of cores detected because on CPUs that support multi- or hyperthreading
  the number of cores reported by \code{detectCores()} is double the number of physical cores. Setting \code{nthreads = 1} will turn off threading. Different C++ routines are used for the matrix fill in the unthreaded (\code{\link{zpm}}) and threaded (\code{\link{zpmCP}}) cases. Speed improvements, if any, may vary.

  \code{sigma_ed} and \code{qt_ed} control the smoothing and threshold parameters passed to the automated edge detection function \code{\link{circle.auto}}

  If you don't like the default color palette there are many other choices. If you like rainbows \code{\link{rygcb}}
  defined in this package produces a relatively perceptually uniform version that's well suited for display on an RGB monitor.
  
}
\value{
  A named list with the current values of the arguments.
}
\author{M.L. Peck \email{mlpeck54@gmail.com}}
\keyword{list}
