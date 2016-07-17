\name{pcafit}
\alias{pcafit}
\title{Principal Component Analysis for Phase Shifting Interferometry}
\description{
  Implements an algorithm proposed by Vargas et al. (2011) to
  estimate phase in PSI with principal components.
}
\usage{
pcafit(images, cp=NULL,
	BGSUB = TRUE, diagpos = c("v", "u", "g"), trace=1, REFINE = TRUE, 
        fringescale=1, zlist=makezlist(), zc0=c(1:3, 6:7), 
	satarget=c(0,0), astig.bath=c(0,0), 
	puw.alg="qual", uselm=FALSE, sgs=1, plots=TRUE, CROP=FALSE)
}
\arguments{
  \item{images}{An array containing the interferogram images}
  \item{cp}{A list describing the pupil boundary, as returned by \code{\link{pupil.pars}}}
  \item{BGSUB}{Logical: subtract pixel-wise mean?}
  \item{diagpos}{Algorithm variant -- see details}
  \item{trace}{Output summary every trace iterations}
  \item{REFINE}{Logical: calculate phase again after finding interferogram edge?}
  \item{fringescale}{Fringe spacing, in waves. Use 1 for single pass, 0.5 for double pass, etc.}
  \item{zlist}{Indexes of Zernike polynomials to fit to wavefront}
  \item{zc0}{Indexes of Zernike coefficients to be removed from net wavefront}
  \item{satarget}{Target 4th and 6th order SA coefficients in non-null tests of aspheres}
  \item{astig.bath}{Astigmatism coefficients for Bath geometry}
  \item{puw.alg}{Algorithm to use for phase unwrapping}
  \item{uselm}{Boolean: use lm() for least squares fit}
  \item{sgs}{Sample Grid Spacing for least squares fits to wavefront values}
  \item{plots}{Logical: plot progress?}
  \item{CROP}{Center and crop maps?}
}
\details{
  \code{images} is a 3 dimensional array with dimensions \code{nrow x ncol x length(phases)},
  where \code{nrow} and \code{ncol} are the number of rows and columns in the individual
  interferogram images.
  
  If \code{cp == NULL} \code{\link{circle.pars}} is called to construct the
  interferogram mask automatically.
  
  \code{diagpos} controls the operation of the underlying algorithm. If \code{diagpos = "u"}
  singular values are ``grouped'' with the matrix \bold{U} in the singular value decomposition
  of the image data. This is the algorithm as originally described in Vargas et al. \code{diagpos = "v"}
  groups singular values with the matrix \bold{V} in the SVD. This is the default for
  consistency with the original implementation by the author. Finally, \code{diagpos = "g"} calls
  the generalized PCA algorithm \code{\link{gpcapsi}} developed by the author.

  If \code{REFINE == TRUE} \code{pcapsi} is called again with the masked image
  data after the calculation of \code{cp}. If \code{!is.null(cp)} in the function
  call \code{REFINE} is ignored.
  
  See \code{\link{wf.net}} for details of the process of creating net and
  smoothed wavefronts from raw unwrapped wavefront maps.
  
  \code{puw.alg} specifies the unwrapping algorithm. You can specify an algorithm
  by choosing \code{puw.alg=c("brcut", "lp", "qual")}.
  
}
\value{
  A list with the following components:
  \item{phase}{Raw, wrapped phase map}
  \item{mod}{The estimated modulation}
  \item{phases}{final phase shift estimates}
  \item{cp}{A list describing the pupil boundary}
  \item{cp.orig}{Uncropped value of cp}
  \item{wf.net}{Net unsmoothed wavefront; a matrix of class \code{"\link{pupil}"}}
  \item{wf.smooth}{Net smoothed wavefront}
  \item{wf.residual}{Difference between net wavefront and polynomial fit}
  \item{fit}{Return value from \code{\link{fitzernikes}}}
  \item{zcoef.net}{Net Zernike coefficients from fit}
  \item{snr}{An estimate of the signal to noise ratio of the raw data}
  \item{eigen}{The singular values of the input crossproduct matrix}
}
\references{
J. Vargas, J. Antonio Quiroga, and T. Belenguer, 
``Phase-shifting interferometry based on principal component analysis,''
\emph{Opt. Lett.} \bold{36}, 1326-1328 (2011)
\href{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-8-1326}{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-8-1326}

J. Vargas, J. Antonio Quiroga, and T. Belenguer, 
``Analysis of the principal component algorithm in phase-shifting interferometry,'' 
\emph{Opt. Lett.} \bold{36}, 2215-2217 (2011)
\href{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-12-2215}{http://www.opticsinfobase.org/ol/abstract.cfm?URI=ol-36-12-2215}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}

\examples{
    fpath <- file.path(find.package(package="zernike"), "psidata")
    files <- scan(file.path(fpath, "files_pca.txt"), what="character")
    for (i in 1:length(files)) files[i] <- file.path(fpath, files[i])
        
        # load the images into an array
        
    images <- load.images(files)

    # parameters for this run
        
    source(file.path(fpath, "parameters.txt"))
        
        # target SA coefficients for numerical null.
        
    sa.t <- 2*zconic(diam,roc,lambda=wavelength*1e-6)
        
      ## default arguments; edge detection works, so use it in subsequent calls
    
    pcfit.v <- pcafit(images, satarget=sa.t)
        
    ## the other algorithm options
    
    pcfit.u <- pcafit(images, cp=pcfit.v$cp, diagpos="u", satarget=sa.t)
    pcfit.g <- pcafit(images, cp=pcfit.v$cp, diagpos="g", satarget=sa.t)
    
    plotn(pcfit.v, pcfit.u, pcfit.g)
    plotn(pcfit.v, pcfit.u, pcfit.g, wftype="smooth")
    
      ## clean up
    
    rm(fpath, files, images, sa.t, diam, roc, wavelength, ps.dir, frames.per.cycle,
       beam.sep, i, pcfit.v, pcfit.u, pcfit.g)
}

\seealso{
  \code{\link{psifit}},
  \code{\link{itfit}},
  \code{\link{wf.net}},
  \code{\link{circle.pars}}.
}
\keyword{models}
