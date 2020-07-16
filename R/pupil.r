## create a unit aperture inside a matrix of arbitrary size
## note the parameter cp is a list with components (xc, yc, rx, ry, obstruct) as returned by pupil.pars()

pupil <- function(zcoef=NULL, zlist=makezlist(), phi=0, piston=0,
                  nrow=512, ncol=nrow, cp=list(xc=256,yc=256,rx=255,ry=255,obstruct=0)) {
  
  prt <- pupil.rhotheta(nrow, ncol, cp)
  rho <- prt$rho
  theta <- prt$theta
  rho.v <- rho[!is.na(rho)]
  theta.v <- theta[!is.na(rho)]
  wf.v <- numeric(length(rho.v))
  if (!is.null(zcoef)) {
    wf.v <- zpm(rho.v, theta.v, phi, maxorder= max(zlist$n)) %*% 
    c(0, zcoef)
  }
  wf <- matrix(NA, nrow=nrow, ncol=ncol)
  wf[!is.na(rho)] <- wf.v
  wf <- wf+piston
  class(wf) <- "pupil"
  wf
}
                  
pupil.arb <- function(zcoef=NULL, zlist=makezlist(), phi=0, piston=0,
	nrow=512, ncol=nrow, cp=list(xc=256,yc=256,rx=255,ry=255,obstruct=0)) {

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
  class(wf) <- "pupil"
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
