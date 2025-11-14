########
## classical fft fringe analysis
########


.up2 <- function(nr, nc=nr) 2^(ceiling(log2(max(nr,nc))))

## FFT fit routine


fftfit <- function(img, cp=NULL,
                   sl=c(1,1), filter=NULL, taper=2, options = zernike_options()) {
  
  if (!is.null(options$use_fftw) && options$use_fftw) {
    fft <- zernike::fft_fftw
    ifft <- zernike::ifft_fftw
  }
  
  nr <- nrow(img)
  nc <- ncol(img)
  npad <- nextn(max(nr,nc))
  if (!is.null(cp)) {
    prt <- pupil.rhotheta(nr, nc, cp)
    img[is.na(prt$rho)] <- mean(img)
  }
  img <- img-mean(img)
  im <- padmatrix(img, npad=npad, fill=0)
  im.fft <- fftshift(fft(im))
  if (is.null(filter)) {
    sldata <- pick.sidelobe(img)
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
  im.fft <- im.fft*gblur(sfilter, sigma=taper)
  sl.fft <- matrix(0, npad, npad)
  xmin <- max(1-sl[1], 1)
  xmax <- min(npad-sl[1], npad)
  ymin <- max(1-sl[2], 1)
  ymax <- min(npad-sl[2], npad)
  sl.fft[xmin:xmax, ymin:ymax] <- im.fft[(sl[1]+(xmin:xmax)),(sl[2]+(ymin:ymax))]
  if (options$plots) {
    plot.cmat(submatrix(sl.fft, size=npad/2))
  }
  cphi <- ifft(fftshift(sl.fft))[1:nr, 1:nc]
  phi <- Arg(cphi)
  mod <- Mod(cphi)
  mod <- mod/max(mod)
  wfnets <- wf_net(phi, mod, cp, options)
  outs0 <- list(rundate = date(), algorithm="classical FFT")
  outs <- c(outs0, wfnets)
  class(outs) <- c(class(outs), "wf_zfit")
  outs
}
                   
                   
