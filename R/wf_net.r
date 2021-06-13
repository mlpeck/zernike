## calculate net and zernike fit wavefronts from any of psifit, fftfit, vortexfit

wf_net <- function(wf.raw, cp, options) {
  zlist <- makezlist(maxorder=options$maxorder)
  nr <- nrow(wf.raw)
  nc <- ncol(wf.raw)
  prt <- pupil.rhotheta(nr, nc, cp)
  rho <- prt$rho
  theta <- prt$theta
  if (options$sgs > 1) {
    subs <- matrix(FALSE, nr, nc)
    subs[seq(1, nr, by=options$sgs), seq(1, nc, by=options$sgs)] <- TRUE
    subs[is.na(rho)] <- FALSE
  } else {
    subs <- !is.na(rho)
  }
  wf.v <- wf.raw[subs]
  rho.v <- rho[subs]
  th.v <- theta[subs]
  rho.v <- rho.v[!is.na(wf.v)]
  th.v <- th.v[!is.na(wf.v)]
  wf.v <- wf.v[!is.na(wf.v)]
  fit <- fitzernikes(wf.v, rho.v, th.v, maxorder=options$maxorder, nthreads=options$nthreads, uselm=options$uselm)
  if (options$uselm) {
    cfit <- coef(fit)
  } else {
    cfit <- fit
  }
  if (sign(cfit[9])*sign(options$satarget[1]) < 0) {
    cfit <- -cfit
    wf.raw <- -wf.raw
    if (!options$uselm) fit <- -fit
  }
  zc.low <- rep(0,15)
  zc.low[options$zc0] <- cfit[options$zc0+1]
  zc.low[c(8,15)] <- zc.low[c(8,15)] + options$satarget
  zc.low[4:5] <- zc.low[4:5] + options$astig.bath
  wf.net <- wf.raw - pupil(zcoef=zc.low, zlist=makezlist(2,6), piston=cfit[1], 
                           nrow=nr, ncol=nc, cp=cp)
  if (options$plots) {
    if (tolower(.Platform$OS.type) == "windows") {
      windows(width=18, height=6)
    } else {
        x11(width=18, height=6)
    }
    split.screen(figs=c(1,3))
    screen(1)
    plot(wf.net, cp=cp, col=options$colors, addContours=FALSE)
    mtext(paste("RMS = ", format(pupilrms(wf.net),digits=3)))
  }
  zcoef.net <- cfit[-1]
  zcoef.net[1:15] <- zcoef.net[1:15] - zc.low
  wf.smooth <- pupil(zcoef=zcoef.net, zlist=zlist, cp=cp, nrow=nr, ncol=nc)
  if (options$plots) {
    screen(2)
    plot(wf.smooth, cp=cp, col=options$colors)
    mtext(paste("RMS = ", format(pupilrms(wf.smooth),digits=3)))
  }
  wf.residual <- wf.net - wf.smooth
  if (options$plots) {
    screen(3)
    plot(wf.residual, cp=cp, col=grey256, addContours=FALSE)
    mtext(paste("RMS = ", format(pupilrms(wf.residual),digits=3)))
    close.screen(all.screens=TRUE)
  }
  list(wf.net=wf.net, wf.smooth=wf.smooth, wf.residual=wf.residual, 
       fit=fit, zcoef.net=zcoef.net)
}

#' Methods for class "wf_fitted"
#' 
#' Summary, print, and plot methods for the returned list of values
#'  from [psifit()], [fftfit()], or [vortexfit()]
#'
#' @param wffit the return values from one of the fringe analysis routines
#' @para digits number of digits to display in print and summary methods
#' @param ... values passed to [plot.pupil()]
#' @return print method returns data frame with Zernike coefficients
summary.wf_fitted <- function(wffit, digits=3) {
  cat("Image size(s)     : ", nrow(wffit$wf.smooth) "x", ncol(wffit$wf.smooth), "\n")
  cat("Unsmoothed RMS    : ", format(pupilrms(wffit$wf.net), digits=digits), "\n")
  cat("Zernike fit RMS   : ", format(sqrt(crossprod(wffit$zcoef.net)), digits=digits), "\n")
  cat("Zernike fit Strehl: ", format(strehlratio(sqrt(crossprod(wffit$zcoef.net))), digits=digits), "\n")
  cat("Zernie fit P-V    : ", format(pupilpv(wffit$wf.smooth), digits=digits), "\n")
  cat("PVr               : ", format(PVr(wffit$wf.smooth, wffit$wf.residual), digits=digits), "\n")
}

print.wf_fitted <- function(wffit, digits=3) {
  if (class(wffit$fit) == "lm") {
    fit <- coef(wffit$fit)
  } else {
    fit <- coef(wffit$fit)
  }
  nz <- length(wffit$zcoef.net)
  znames <- paste("Z", 0:nz, sep="")
  df.zernikes <- data.frame(Z = znames, Zcoef.Raw = fit, Zcoef.net = c(0, wffit$zcoef.net))
  summary.wf_fitted(wffit, digits=digits)
  cat("\n")
  print(df.zernikes, digits=digits, row.names=FALSE)
  df.zernikes
}

plot.wf_fitted <- function(wffit, wftype="smooth", ...) {
  wf <- get(paste(wf, wftype, sep="."), wffit)
  plot.pupil(wf, cp=wffit$cp, ...)
}
  
