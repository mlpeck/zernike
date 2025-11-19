########
## classical fft fringe analysis
########


.up2 <- function(nr, nc=nr) 2^(ceiling(log2(max(nr,nc))))

## FFT fit routine


fftfit <- function(img, cp=NULL,
                   sl=NULL, filter=NULL,
                   options = zernike_options()) {
  
  if (!is.null(options$use_fftw) && options$use_fftw) {
    fft <- zernike::fft_fftw
    ifft <- zernike::ifft_fftw
  }
  
  nr <- nrow(img)
  nc <- ncol(img)

  if (!is.null(cp)) {
    prt <- pupil.rhotheta(nr, nc, cp)
    img[is.na(prt$rho)] <- mean(img)
  }
  img <- img-mean(img)

  if (nextn(nr) != nr || nextn(nc) != nc) {
    npad <- nextn(max(nr,nc))
    nr <- nc <- npad
    img <- padmatrix(img, npad=npad)
  }

  if (is.null(filter) || is.null(sl)) {
    sldata <- pick.sidelobe(img)
    sl <- sldata$sl
    filter <- sldata$filter
  }
  xs <- (1:nr) - (nr %/% 2 + 1)
  ys <- (1:nc) - (nc %/% 2 + 1)
  carrier <- outer(sl[1]*xs/nr, sl[2]*ys/nc, function(x,y) exp(-2*pi*1i*(x+y)))
  im.fft <- fftshift(fft_cx(carrier*img))

  step <- outer(xs, ys, function(x,y) as.numeric(x^2+y^2 <= filter^2))
  im.fft <- im.fft*step

  cphi <- ifft(ifftshift(im.fft))[1:nrow(img), 1:ncol(img)]
  phi <- Arg(cphi)
  mod <- Mod(cphi)
  mod <- mod/max(mod)
  wfnets <- wf_net(phi, mod, cp, options)
  outs0 <- list(rundate = date(), algorithm="classical FFT", sl=sl, filter=filter)
  outs <- c(outs0, wfnets)
  class(outs) <- c(class(outs), "wf_zfit")
  outs
}
                   
                   
