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
As of ver. 3.7.2 the values of \code{options$isoseq} and \code{cp$obstruct} determine what is called to fit Zernikes.
If \code{isoseq} is \code{TRUE} it's either \code{\link{zpm_cart}} or \code{\link{zapm_cart}}, with the latter
called automatically if \code{cp$obstruct > 0}. Otherwise the calls are made to
\code{\link{zpm}} or \code{\link{zapmC}}.
}

\examples{

## illustrates effect of incremental addition of an obstruction
## and standalone use of wf_net
  
require(zernike)
fpath <- file.path(find.package(package="zernike"), "psidata")
files <- scan(file.path(fpath, "files.txt"), what="character")
for (i in 1:length(files)) files[i] <- file.path(fpath, files[i])

# load the images into an array

images <- load.images(files)

# parameters for this run

source(file.path(fpath, "parameters.txt"))

# phase shifts

phases <- wrap((0:(dim(images)[3]-1))/frames.per.cycle*2*pi)
phases <- switch(ps.dir, ccw = -phases, cw = phases, phases)

# target SA coefficients for numerical null.

sa.t <- sconic(diam,roc,lambda=wavelength)
zopt <- psfit_options()
zopt$satarget <- sa.t
psfit <- psifit(images, phases, psialg="ls", options=zopt)

# get back the raw wavefront

wf.raw <- qpuw(psfit$phi, psfit$mod)
cp <- psfit$cp

eps <- seq(0.1, 0.5, by=0.1)
ne <- length(eps)
pad0 <- rep(0, ne)

## collect some basic data

df.annfits <- data.frame(eps=c(0, eps), sa.null = c(sa.t[1], pad0), sa.obs = c(psfit$fit[9], pad0), sa.net = c(psfit$zcoef.net[8], pad0),
                         rms.net = c(sqrt(crossprod(psfit$zcoef.net)), pad0), defocus=c(psfit$fit[4], pad0))
for (i in seq_along(eps)) {
  cp$obstruct <- eps[i]
  sa.t <- sconic(diam, roc, eps=eps[i], lambda=wavelength)
  zopt$satarget <- sa.t
  fiti <- wf_net(wf.raw, cp=cp, options=zopt)
  df.annfits[i+1, 2:ncol(df.annfits)] <- c(sa.t[1], fiti$fit[9], fiti$zcoef.net[8], sqrt(crossprod(fiti$zcoef.net)), fiti$fit[4])
}

rm(eps)
attach(df.annfits)

X11()
plot(eps, sa.null, type="b", ylim=range(df.annfits[,2:ncol(df.annfits)]), xlim=c(0, 0.8), main="Trends with obstruction fraction")
points(eps, sa.obs, type="b", pch=2, col=2, lty=2)
points(eps, sa.net, type="b", pch=3, col=3, lty=3)
points(eps, rms.net, type="b", pch=4, col=4, lty=4)
points(eps, defocus, type="b", pch=5, col=5, lty=5)
grid()
legend(x=0.6, y=0, legend=names(df.annfits)[-1], col=1:5, lty=1:5, pch=1:5)

## plot the first and last smoothed wavefronts to show difference is completely symmetrical

plotn(psfit, fiti, wftype="smooth", qt=c(0,1))
detach(df.annfits)
print(df.annfits, digits=2)

}
\keyword{Graphics}% use one of  RShowDoc("KEYWORDS")
\keyword{smooth}% __ONLY ONE__ keyword per line
