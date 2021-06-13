########
## Utilities for fft fringe analysis
########


.up2 <- function(nr, nc=nr) 2^(ceiling(log2(max(nr,nc))))

## FFT fit routine
#
## parameters:
## FFT fit routine
#
## parameters:
#
#	imagedata: matrix containing greyscale values of interferogram image
#	cp: a list with components (xc, yc, rx, ry, obstruct) describing pupil parameters
#	fringescale: scale factor for fringes: 1 for single pass, .5 for double
#	sl: approximate location of desired sidelobe in the form c(x,y)
#	filter: size of filter around DC
#	taper: amount to taper edge of half plane cut
#	zlist: Zernikes to fit. Defaults to order 14
#	zc0: Zernikes we want zeroed. Defaults to tilts, defocus, and coma
#	satarget: SA term for desired asphere, if appropriate
#	astig.bath: astigmatism due to Bath geometry
#	puw.alg: phase unwrapping algorithm



fftfit <- function(imagedata, cp=NULL, 
                   sl=c(1,1), filter=NULL, taper=2, options = psfit_options()) {
  
  nr <- nrow(imagedata)
  nc <- ncol(imagedata)
  npad <- nextn(max(nr,nc))
  if (!is.null(cp)) {
    prt <- pupil.rhotheta(nr, nc, cp)
    imagedata[is.na(prt$rho)] <- mean(imagedata)
  }
  imagedata <- imagedata-mean(imagedata)
  im <- padmatrix(imagedata, npad=npad, fill=0)
  im.fft <- fftshift(fft(im))
  if (is.null(filter)) {
    sldata <- pick.sidelobe(imagedata)
    sl <- sldata$sl
    filter <- sldata$filter
  }
  if (filter > 0) {
    xs <- 2*(-(npad/2):(npad/2-1))/filter
    rho2 <- outer(xs, xs, function(x,y) x^2+y^2)
    im.fft <- im.fft*(1-exp(-rho2/2))
  }
  sl <- round(sl) #round sidelobe position to integer values
  ## decide which half plane to 0
  sfilter <- matrix(1,npad,npad)
  if (sl[1] == 0) {
    if (sl[2] > 0) {
      sfilter[, 1:(npad/2+1)] <- 0
    } else {
      sfilter[, (npad/2+1):npad] <- 0
    }
  } else if (sl[2] == 0) {
    if (sl[1] > 0) {
      sfilter[1:(npad/2+1),] <- 0
    } else {
      sfilter[(npad/2+1):npad,] <- 0
    }
  } else {
    dydx <- -sl[1]/sl[2]
    xcut <- (1:npad)-(npad/2+1)
    ycut <- (npad/2+1) + round(dydx*xcut)
    ycut[ycut<1] <- 1
    ycut[ycut>npad] <- npad
    if (sl[2] > 0) {
      for (i in 1:npad) sfilter[i, 1:ycut[i]] <- 0
    } else {
      for (i in 1:npad) sfilter[i, ycut[i]:npad] <- 0
    }
  }
  im.fft <- im.fft*gblur(sfilter, fw=taper)
  sl.fft <- matrix(0, npad, npad)
  xmin <- max(1-sl[1], 1)
  xmax <- min(npad-sl[1], npad)
  ymin <- max(1-sl[2], 1)
  ymax <- min(npad-sl[2], npad)
  sl.fft[xmin:xmax, ymin:ymax] <- im.fft[(sl[1]+(xmin:xmax)),(sl[2]+(ymin:ymax))]
  if (options$plots) {
    plot.cmat(submatrix(sl.fft, size=npad/2))
  }
  cphase <- fft(fftshift(sl.fft), inv=TRUE)[1:nr, 1:nc]
  phase <- Arg(cphase)
  mod <- Mod(cphase)
  mod <- mod/max(mod)
  if (is.null(cp)) {
    cp <- circle.pars(mod, plot=options$plots)
  }
  cp.orig <- cp
  if (options$crop) {
    mod <- crop(mod, cp)$im
    phase <- crop(phase, cp)
    cp <- phase$cp
    phase <- phase$im
  }
  nr <- nrow(phase)
  nc <- ncol(phase)
  prt <- pupil.rhotheta(nr, nc, cp)
  phase[is.na(prt$rho)] <- NA
  class(phase) <- "pupil"
  wf.raw <- switch(options$puw_alg,
                   qual = qpuw(phase, mod),
                   brcut = brcutpuw(phase),
                   qpuw(phase, mod)
  )
  wf.raw <- options$fringescale*wf.raw
  class(wf.raw) <- "pupil"
  wf.nets <- wf_net(wf.raw, cp, options)
  outs <- list(phi=phase, mod=mod, cp=cp,cp.orig=cp.orig,
       wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth, 
       wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net)
  class(outs) <- "wf_fitted"
  outs
}
                   
                   
