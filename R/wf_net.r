## calculate net and zernike fit wavefronts from any of psifit, fftfit, vortexfit

wf_net <- function(wf.raw, cp, options) {

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
  fit <- fitzernikes(wf.v, rho.v, th.v, 
                     maxorder=options$maxorder, 
                     nthreads=options$nthreads, 
                     uselm=options$uselm,
                     isoseq=options$isoseq)
  if (options$uselm) {
    cfit <- coef(fit)
  } else {
    cfit <- fit
  }
  if (isoseq) {
    ind.ptf <- c(1:3, 5)
    ind.sa4 <- 13
    ind.sa6 <- 25
    ind.astig <- c(6, 4)
    ind.coma <- 8:9
    zc.low <- numeric(28)
  } else {
    ind.ptf <- 1:4
    ind.sa4 <- 9
    ind.sa6 <- 16
    ind.astig <- 5:6
    ind.coma <- 7:8
    zc.low <- numeric(16)
  }
  if (sign(cfit[ind.sa4])*sign(options$satarget[1]) < 0) {
    cfit <- -cfit
    wf.raw <- -wf.raw
    if (!options$uselm) fit <- -fit
  }
  
  zc.low[ind.ptf] <- cfit[ind.ptf]
  zc.low[c(ind.sa4, ind.sa6)] <- zc.low[c(ind.sa4, ind.sa6)] + options$satarget
  zc.low[ind.astig] <- zc.low[ind.astig] + options$astig.bath
  if (is.element(6, options$zc0)) {
    zc.low[ind.coma] <- cfit[ind.ptf]
  }
  wf.net <- wf.raw - pupil(zcoef=zc.low, maxorder=6, isoseq=options$isoseq,
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
  zcoef.net[1:length(zc.low)] <- zcoef.net[1:length(zc.low)] - zc.low
  wf.smooth <- pupil(zcoef=zcoef.net, maxorder=options$maxorder, isoseq=options$isoseq, 
                     cp=cp, nrow=nr, ncol=nc)
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
#' @param digits number of digits to display in print and summary methods
#' @param ... values passed to [plot.pupil()]
#' @return print method returns data frame with Zernike coefficients
summary.wf_fitted <- function(wffit, digits=3) {
  cat("Image size(s)     : ", nrow(wffit$wf.smooth), " x ", ncol(wffit$wf.smooth), "\n")
  cat("Unsmoothed RMS    : ", format(pupilrms(wffit$wf.net), digits=digits), "\n")
  cat("Zernike fit RMS   : ", format(sqrt(crossprod(wffit$zcoef.net)), digits=digits), "\n")
  cat("Zernike fit Strehl: ", format(strehlratio(sqrt(crossprod(wffit$zcoef.net))), digits=digits), "\n")
  cat("Zernie fit P-V    : ", format(pupilpv(wffit$wf.smooth), digits=digits), "\n")
  cat("PVr               : ", format(PVr(wffit$wf.smooth, wffit$wf.residual), digits=digits), "\n")
}

print.wf_fitted <- function(wffit, digits=3) {
  if (is.element("lm", class(wffit$fit))) {
    fit <- coef(wffit$fit)
  } else {
    fit <- wffit$fit
  }
  nz <- length(wffit$zcoef.net)
  znames <- paste("Z", 0:nz, sep="")
  df.zernikes <- data.frame(Z = znames, Zcoef.Raw = fit, Zcoef.net = c(0, wffit$zcoef.net))
  summary.wf_fitted(wffit, digits=digits)
  cat("\n")
  print(df.zernikes, digits=digits, row.names=FALSE)
  invisible(df.zernikes)
}

plot.wf_fitted <- function(wffit, wftype="smooth", ...) {
  wf <- get(paste("wf", wftype, sep="."), wffit)
  plot.pupil(wf, cp=wffit$cp, ...)
}
  
