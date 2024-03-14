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
                     eps=cp$obstruct,
                     maxorder=options$maxorder, 
                     nthreads=options$nthreads, 
                     uselm=options$uselm,
                     isoseq=options$isoseq,
                     usecirc=options$usecirc)
  if (options$uselm) {
    cfit <- coef(fit)
  } else {
    cfit <- fit
  }
  if (options$isoseq) {
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
    zc.low[ind.coma] <- cfit[ind.coma]
  }
  wf.net <- wf.raw - pupil(zcoef=zc.low, maxorder=6, 
                           isoseq=options$isoseq, usecirc=options$usecirc,
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
  zcoef.net <- cfit
  zcoef.net[1:length(zc.low)] <- zcoef.net[1:length(zc.low)] - zc.low
  wf.smooth <- pupil(zcoef=zcoef.net, maxorder=options$maxorder, 
                     isoseq=options$isoseq, usecirc=options$usecirc,
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
  outs <- list(wf.net=wf.net, wf.smooth=wf.smooth, wf.residual=wf.residual, 
       fit=fit, zcoef.net=zcoef.net[-1])
  class(outs) <- c(class(outs), "wf_zfit")
  outs
}

#' Methods for class "wf_zfit"
#' 
#' Summary, print, plot, and invert methods for the returned list of values
#'  from [psifit()], [fftfit()], or [vortexfit()]
#'
#' @param wffit the return values from one of the fringe analysis routines or [wf_net()]
#' @param digits number of digits to display in print and summary methods
#' @param printnow send output to console?
#' @param ... values passed to [plot.pupil()]
#'
#' @details The invert method negates the values of wavefronts and Zernike coefficients and returns the adjusted input
#'
#' @return summary and print methods return data frame with wavefront summaries and Zernike coefficients
summary.wf_zfit <- function(wffit, digits=3, printnow=TRUE) {
  df.sum <- data.frame(Name    = c("Current time      : ",
                                   "Run time          : ",
                                   "Algorithm         : ",
                                   "Image size(s)     : ",
                                   "Unsmoothed RMS    : ",
                                   "Zernike fit RMS   : ",
                                   "Zernike fit Strehl: ",
                                   "Zernike fit P-V   : ",
                                   "PVr               : "),
                       Value  = c(date(),
                                  if (!is.null(wffit$rundate)) wffit$rundate else "",
                                  if (!is.null(wffit$algorithm)) wffit$algorithm else "unknown",
                                  paste(nrow(wffit$wf.smooth), "x", ncol(wffit$wf.smooth)),
                                  format(pupilrms(wffit$wf.net),digits=digits),
                                  format(hypot(wffit$zcoef.net),digits=digits),
                                  format(strehlratio(hypot(wffit$zcoef.net)),digits=digits),
                                  format(pupilpv(wffit$wf.smooth),digits=digits),
                                  format(PVr(wffit$wf.smooth, wffit$wf.residual),digits=digits))
  )
  if (printnow) {
    print(df.sum, digits=digits, right=FALSE, justify="left", row.names=FALSE)
  }
  invisible(df.sum)
}

print.wf_zfit <- function(wffit, digits=3, abnames=TRUE, printnow=TRUE) {
  if (is.element("lm", class(wffit$fit))) {
    fit <- coef(wffit$fit)
  } else {
    fit <- wffit$fit
  }
  nz <- length(fit)
  sqnz <- sqrt(nz)
  if (sqnz %% 1 > sqrt(.Machine$double.eps)) {
    ## Zernikes are (probably) in ISO sequence
    maxorder <- (-3 + sqrt(1+8*nz))/2
    zlist <- makezlist.iso(maxorder=maxorder)
  } else {
    maxorder <- 2*(sqnz-1)
    zlist <- makezlist(maxorder=maxorder)
  }
  znames <- paste("Z", 0:(nz-1), sep="")
  mmult <- 1 - 2*as.numeric(zlist$t == "s")
  n <- zlist$n
  m <- zlist$m * mmult
  if (abnames) {
    ab <- rep("", nz)
    ab[n==0 & m==0] <- "piston"
    ab[n==2 & m==0] <- "defocus"
    ab[n>=4 & m==0] <- paste("spherical ", n[n>=4 & m==0], "'th", sep="")
    
    ab[n==1 & m==1] <- "x tilt"
    ab[n==1 & m==-1] <- "y tilt"
    cond <- n>=3 & m==1
    ab[cond] <- paste("x coma ", n[cond] + 1, "'th", sep="")
    cond <- n>=3 & m== -1
    ab[cond] <- paste("y coma ", n[cond] + 1, "'th", sep="")
    
    ab[m==2] <- paste("cos astigmatism ", n[m==2] + 2, "'th", sep="")
    ab[m== -2] <- paste("sin astigmatism ", n[m== -2] + 2, "'th", sep="")
    
    ab[m==3] <- paste("cos trefoil ", n[m==3] + 3, "'th", sep="")
    ab[m== -3] <- paste("sin trefoil ", n[m== -3] + 3, "'th", sep="")
    
    df.zernikes <- data.frame(Z = znames, n=n, m=m, aberration=ab,
                              zcoef.raw = fit, zcoef.net = c(0, wffit$zcoef.net))
  } else {
    df.zernikes <- data.frame(Z = znames, n=n, m=m,
                              zcoef.raw = fit, zcoef.net = c(0, wffit$zcoef.net))
  }
  if (printnow) {
    summary(wffit, digits=digits)
    cat("\n")
    print(df.zernikes, digits=digits, row.names=FALSE)
  }
  invisible(df.zernikes)
}

plot.wf_zfit <- function(wffit, wftype="smooth", ...) {
  wf <- get(paste("wf", wftype, sep="."), wffit)
  plot.pupil(wf, cp=wffit$cp, ...)
}

invert <- function(wffit) UseMethod("invert", wffit)

invert.wf_zfit <- function(wffit) {
  if (exists("phi", wffit)) {
    wffit$phi <- -wffit$phi
  }
  if (exists("phases", wffit)) {
    wffit$phases <- -wffit$phases
  }
  wffit$wf.net <- -wffit$wf.net
  wffit$wf.smooth <- -wffit$wf.smooth
  wffit$wf.residual <- -wffit$wf.residual
  wffit$zcoef.net <- -wffit$zcoef.net
  if (!is.element("lm", class(wffit$fit))) {
    wffit$fit <- -wffit$fit
  }
  wffit
}

