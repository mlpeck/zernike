\name{psifit}
\alias{psifit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phase Shifting Interferometry}
\description{
  High level function for Least squares analysis of
  phase shifted interferograms.
}
\usage{
psifit(images, phases, cp = NULL, satarget = NULL, psialg = "ls", options = psfit_options)
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
  
  \item{ls}{least squares with known phase shifts}
  \item{aia}{the ``advanced iterative algorithm`` \code{\link{aiapsi}}}
  \item{pc1}{pca with group_diag = "v"}
  \item{pc2}{pca with group_diag = "u"}
  \item{gpc}{my generalized PC algorithm in \code{\link{gpcapsi}}}
  \item{gpcthentilt}{first \code{\link{gpcapsi}} the \code{\link{tiltpsi}}}
  \item{tilt}{tiltpsi}
  
}
\value{
A list with the following components
\item{phi}{wrapped phase estimate}
\item{mod}{modulation estimate}
\item{phases}{phase shifts}
\item{cp}{the interferogram boundary}
\item{wf.net}{net, unfiltered wavefront (see \code{\link{wf_net}})
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
\keyword{Models}
