## create a unit aperture inside a matrix of arbitrary size
## note the parameter cp is a list with components (xc, yc, rx, ry, obstruct) as returned by pupil.pars()

pupil <- function(zcoef=NULL, maxorder=12L, isoseq=FALSE, 
                  phi=0, piston=NULL,
                  nrow=640, ncol=nrow, 
                  cp=list(xc=320.5,yc=320.5,rx=319.5,ry=319.5,obstruct=0)) {
  
  prt <- pupil.rhotheta(nrow, ncol, cp)
  rho <- prt$rho
  theta <- prt$theta
  rho.v <- rho[!is.na(rho)]
  theta.v <- theta[!is.na(rho)] - pi*phi/180
  wf.v <- numeric(length(rho.v))
  if (!is.null(zcoef)) {
    if (!is.null(piston)) {
      zcoef <- c(piston, zcoef)
    }
    if (isoseq) {
      wf.v <- zpm_cart(x=rho.v*cos(theta.v), y=rho.v*sin(theta.v), maxorder=maxorder) %*% zcoef
    } else {
      wf.v <- zpm(rho.v, theta.v, maxorder = maxorder) %*% zcoef
    }
  }
  wf <- matrix(NA, nrow=nrow, ncol=ncol)
  wf[!is.na(rho)] <- wf.v
  class(wf) <- c("pupil", class(wf))
  wf
}
                  
pupil.arb <- function(zcoef=NULL, zlist=makezlist(), 
                      phi=0, piston=0,
                      nrow=640, ncol=nrow, 
                      cp=list(xc=320.5,yc=320.5,rx=319.5,ry=319.5,obstruct=0)) {

  prt <- pupil.rhotheta(nrow, ncol, cp)
  rho <- prt$rho
  theta <- prt$theta
  rho.v <- rho[!is.na(rho)]
  theta.v <- theta[!is.na(rho)]
  wf.v <- numeric(length(rho.v))
  if (!is.null(zcoef)) {
    wf.v <- zpm.arb(rho.v, theta.v, phi, zlist=zlist) %*% zcoef
  }
  wf <- matrix(NA, nrow=nrow, ncol=ncol)
  wf[!is.na(rho)] <- wf.v
  wf <- wf+piston
  class(wf) <- c("pupil", class(wf))
  wf
}

## estimate of rms over pupil

pupilrms <- function(pupil) {
    sd(as.vector(pupil), na.rm=TRUE)
}

## estimate of p-v over pupil

pupilpv <- function(pupil) {
    max(pupil, na.rm=TRUE) - min(pupil, na.rm=TRUE)
}

#' Zygo's "robust" PV
#'
#' A peak to valley error estimate that reduces the effect of noise and artifacts
#'
#' @param wf.zfit matrix containing the smoothed Zernike fit wavefront
#' @param wf.residual matrix of the difference between the raw wavefront
#'  and the Zernike fit. These values are returned by [wf_net()]
#' @return the estimated PVr
#' @references
#'  Evans, C. (2009) Optical Engineering 48(4), 43605.
#'  <https://doi.org/10.1117/1.3119307>
#'
#' @details
#'  no check is performed on the wavefronts, so it's the user's
#'  responsibility to make sure these come from the same source
PVr <- function(wf.zfit, wf.residual) {
  pvr <- pupilpv(wf.zfit) + 3*sd(wf.residual, na.rm=TRUE)
  pvr
}

## Mahajan's approximation to Strehl ratio

strehlratio <- function(rms) {
    exp(-(2*pi*rms)^2)
}

## a summary method for pupils

summary.pupil <- function(wf) {
	cat("Size:   ",nrow(wf), "x",ncol(wf),"\n")
	cat("RMS    =", format(pupilrms(wf), digits=3), "\n")
	cat("P-V    =", format(pupilpv(wf), digits=3), "\n")
	cat("Strehl =", format(strehlratio(pupilrms(wf)), digits=3), "\n")
}

## polar coordinates of the pupil

pupil.rhotheta <- function(nrow, ncol, cp) {
	xs <- ((1:nrow)-cp$xc)/cp$rx
	ys <- ((1:ncol)-cp$yc)/cp$ry
	rho <- function(x,y) sqrt(x^2+y^2)
	theta <- function(x,y) atan2(y,x)
	rho.mat <- outer(xs, ys, rho)
	theta.mat <- outer(xs,ys,theta)
	rho.mat[rho.mat>1] <- NA
	rho.mat[rho.mat<cp$obstruct] <- NA
	theta.mat[is.na(rho.mat)] <- NA
	list(rho=rho.mat, theta=theta.mat)
}
