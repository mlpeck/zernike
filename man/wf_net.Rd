\name{wf_net}
\alias{wf_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wavefront smoothing}
\description{
  Calculate net and smoothed wavefronts from a raw wavefront containing
  low order nuisance aberrations.
}
\usage{
wf_net(wf.raw, cp, options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wf.raw}{Raw wavefront to be processed}
  \item{cp}{a list describing the pupil boundary}
  \item{options}{a list of options. See \code{\link{psfit_options}}}
}
\details{
Called by \code{\link{psifit}}
}
\value{
  A list with the following components:
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\keyword{Graphics}% use one of  RShowDoc("KEYWORDS")
\keyword{smooth}% __ONLY ONE__ keyword per line
