\name{psifit}
\alias{psifit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phase Shifting Interferometry}
\description{
  High level function for Least squares analysis of
  phase shifted interferograms.
}
\usage{
psifit(images, phases, cp = NULL, satarget = NULL, psialg = "ls", options = psfit_options())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{images}{An array containing the interferogram images}
  \item{phases}{A vector of phase shifts}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{satarget}{Target 4th and 6th order SA coefficients in non-null tests of aspheres}
  \item{psialg}{String identifying the PSI algorithm to use}
  \item{options}{a list of options}
}
\details{
  \code{images} is a 3 dimensional array with dimensions \code{nrow x ncol x length(phases)},
  where \code{nrow} and \code{ncol} are the number of rows and columns in the individual
  interferogram images.
  
  The current values recognized for \code{psialg} are
  
  \describe{
  \item{ls}{least squares with known phase shifts}
  \item{aia}{the ``advanced iterative algorithm`` \code{\link{aiapsi}}}
  \item{pc1}{pca with pcalg="pc1"}
  \item{pc2}{pca with pcalg="pc2"}
  \item{pc3}{pca with pcalg="pc3"}
  \item{gpc}{my generalized PC algorithm in \code{\link{gpcapsi}}}
  \item{pc1thenaia}{first \code{\link{pcapsi}} then \code{\link{aiapsi}}}
  \item{gpcthentilt}{first \code{\link{gpcapsi}} then \code{\link{tiltpsi}}}
  \item{tilt}{tiltpsi}
  }
  
}
\value{
A list with the following components
\item{phi}{wrapped phase estimate}
\item{mod}{modulation estimate}
\item{phases}{phase shifts}
\item{cp}{the interferogram boundary}
\item{wf.net}{net, unfiltered wavefront (see \code{\link{wf_net}})}
\item{wf.smooth}{Zernike fit wavefront}
\item{wf.residual}{the difference}
\item{fit}{Coefficients of Zernike fit to wavefront}
\item{zcoef.net}{Net Zernike coefficients}
\item{extras}{any extra data returned by low level functions}
  
}
\references{
%% ~put references to the literature/web site here ~
}
\author{M.L. Peck \email{mlpeck54@gmail.com}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lspsi}},
  \code{\link{aiapsi}},
  \code{\link{tiltpsi}},
  \code{\link{gpcapsi}},
  \code{\link{pcapsi}}
  
}
\examples{
  ## reuse the files from the demo for an example of two stage fitting
  ## using gpca then tiltpsi
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
zopt$ptol <- 0.01
tfit <- psifit(images, phases, psialg="gpcthentilt", options=zopt)
}

\keyword{Models}
