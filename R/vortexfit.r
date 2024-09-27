##' Vortex transform.
##' 
##' Fringe analysis by Vortex aka Spiral Quadrature transform.
##' 
##' @param imagedata matrix containing the interferogram data
##' @param cp list with circle parameters describing interferogram location. Defaults to NULL
##' @param filter size of filter to remove background
##' @param fw.o size of gaussian blur to smooth orientation estimate
##' @param options A list with general fitting and display options. See [psfit_options()].
##' 
##' @return a list with wavefront estimates, wrapped phase, modulation, etc.
##'         The return has S3 class 'wf_zfit' with plot, print, summary, and report methods.
##' 
##' @seealso This is one of two routines provided for analysis of single interferograms,
##'   along with [fftfit()]. This \emph{may} be suitable for interferograms with
##'   closed fringes.
##'
##' @details Implements the Vortex or spiral phase quadrature transform method
##'  of Larkin et al. (2001) [https://doi.org/10.1364/JOSAA.18.001862] including the 
##'  fringe orientation estimation approach in Larkin (2005) [https://doi.org/10.1364/OPEX.13.008097].
##'  Thanks to Steve Koehler for ideas on implementation details.
##' 
##' @section Warning:
##'   This routine is offered as is with no license, as it may be in violation of one or more
##'   US and international patents.
##'
##' @examples
##' require(zernike)
##' fpath <- file.path(find.package(package="zernike"), "psidata")
##' fname <- "Image197.jpg"
##' img <- load.images(file.path(fpath, fname))
##' 
##' # parameters for this run
##' 
##' source(file.path(fpath, "parameters.txt"))
##' 
##' # target SA coefficients for numerical null.
##' 
##' sa.t <- sconic(diam,roc,lambda=wavelength)
##' zopt <- psfit_options()
##' zopt$satarget <- sa.t
##' 
##' # display an interferogram
##' 
##' if (tolower(.Platform$OS.type) == "windows") windows() else x11()
##' image(1:nrow(img), 1:ncol(img), img, col=grey256, asp=1,
##'  xlab="X", ylab="Y", useRaster=TRUE)
##' mtext("Sample Interferogram")
##' 
##' if (tolower(.Platform$OS.type) == "windows") windows() else x11()
##' vfit <- vortexfit(img, filter=15, fw.o=10, options=zopt)
##'
##' # do "classical FFT" based fit and compare results
##'
##' dev.set(dev.next())
##' ftfit <- fftfit(img, cp=vfit$cp, sl=c(32, 0), filter=15, options=zopt)
##'
##' plotn(ftfit, vfit, labels=c("fft", "vortex"))
vortexfit <- function(imagedata, cp=NULL, filter=NULL, fw.o=10, options=psfit_options()) {
  
  if (require(fftwtools)) {
    fft <- fftwtools::fftw2d
  }
  
  nr <- nrow(imagedata)
  nc <- ncol(imagedata)
  npad <- nextn(max(nr,nc))
  xs <- (-npad/2):(npad/2-1)
  im <- padmatrix(imagedata, npad=npad, fill=mean(imagedata))
  im.fft <- fftshift(fft(im))
  if (is.null(filter)) {
    sldata <- pick.sidelobe(imagedata, logm=TRUE)
    filter <- sldata$filter
  }
  if (filter > 0) {
    xss <- 2*xs/filter
    rho2 <- outer(xss, xss, function(x,y) x^2+y^2)
    im.fft <- im.fft*(1-exp(-rho2/2))
  }
  theta <- function(x,y) atan2(y,x)
  phi <- outer(xs, xs, theta)
  im.nb <- Re(fft(fftshift(im.fft), inverse=TRUE))/(npad^2)
  D <- fft(fftshift(exp(1i*phi)*im.fft), inverse=TRUE)/(npad^2)
  D2 <- fft(fftshift(exp(2*1i*phi)*im.fft), inverse=TRUE)/(npad^2)
  sx <- (D^2-im.nb*D2)[1:nr, 1:nc]
  if (fw.o > 0) {
    rsx <- gblur(Re(sx), fw=fw.o)
    isx <- gblur(Im(sx), fw=fw.o)
    sx <- rsx + 1i*isx
  }
  orient <- Arg(sx)
  mod.o <- Mod(sx)
  if (options$plots) {
    image(1:nr, 1:nc, orient, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE)
    mtext("Orientation")
  }
  dir <- wrap(pi * qpuw(orient, mod.o))
  if (options$plots) {
    image(1:nr, 1:nc, dir, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE)
    mtext("Direction")
  }
  Q <- Re(D[1:nr, 1:nc] * exp(-1i*dir))
  phi <- atan2(Q, im.nb[1:nr, 1:nc])
  mod <- sqrt(Q^2+(im.nb[1:nr, 1:nc])^2)
  mod <- (mod-min(mod))/(max(mod)-min(mod))
  wfnets <- wf_net(phi, mod, cp, options)
  outs0 <- list(rundate = date(), algorithm = "Vortex",
                im.bgclean=im.nb[1:nr, 1:nc], orient=orient, dir=dir)
  outs <- c(outs0, wfnets)
  class(outs) <- c(class(outs), "wf_zfit")
  outs
}
                      
