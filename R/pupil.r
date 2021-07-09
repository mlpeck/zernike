## defaults for matrix size and cp for synthetic wavefront construction

nrow.default <- 640
ncol.default <- nrow.default
cp.default <- list(xc=320.5, yc=320.5, rx=319.5, ry=319.5, obstruct=0)

## create a unit aperture inside a matrix of arbitrary size
## note the parameter cp is a list with components (xc, yc, rx, ry, obstruct) as returned by pupil.pars()

pupil <- function(zcoef=NULL, maxorder=14L, isoseq=FALSE, 
                  phi=0, piston=NULL,
                  nrow=nrow.default, ncol=ncol.default, 
                  cp=cp.default) {
  
  prt <- pupil.rhotheta(nrow, ncol, cp)
  rho <- prt$rho
  theta <- prt$theta
  rho.v <- rho[!is.na(rho)]
  theta.v <- theta[!is.na(rho)] - pi*phi/180
  wf.v <- numeric(length(rho.v))
  if (cp$obstruct > 0) {
    use.circ <- FALSE
  } else {
    use.circ <- TRUE
  }
  if (!is.null(zcoef)) {
    if (!is.null(piston)) {
      zcoef <- c(piston, zcoef)
    }
    if (isoseq) {
      if (use.circ) {
        wf.v <- zpm_cart(x=rho.v*cos(theta.v), y=rho.v*sin(theta.v), maxorder=maxorder) %*% zcoef
      } else {
        wf.v <- zapm_cart(x=rho.v*cos(theta.v), y=rho.v*sin(theta.v), maxorder=maxorder) %*% zcoef
      }
    } else {
      if (use.circ) {
        wf.v <- zpm(rho.v, theta.v, maxorder = maxorder) %*% zcoef
      } else {
        wf.v <- zapmC(rho.v, theta.v, maxorder = maxorder) %*% zcoef
      }
    }
  }
  wf <- matrix(NA, nrow=nrow, ncol=ncol)
  wf[!is.na(rho)] <- wf.v
  class(wf) <- c("pupil", class(wf))
  wf
}
                  
pupil.arb <- function(zcoef=NULL, zlist=makezlist(), 
                      phi=0, piston=0,
                      nrow=nrow.default, ncol=ncol.default, 
                      cp=cp.default) {

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

summary.pupil <- function(wf, digits=3) {
	cat("Size:   ",nrow(wf), "x",ncol(wf),"\n")
	cat("RMS    =", format(pupilrms(wf), digits=digits), "\n")
	cat("P-V    =", format(pupilpv(wf), digits=digits), "\n")
	cat("Strehl =", format(strehlratio(pupilrms(wf)), digits=digits), "\n")
}

## a plot method for pupils

plot.pupil <- function(wf, cp=NULL, col=topo.colors(256), addContours=TRUE, 
                       cscale=FALSE, eqa=FALSE, zlim=NULL, ...) {
    nr <- nrow(wf)
    nc <- ncol(wf)
    if(is.null(zlim)) zlim <- range(wf, finite=TRUE)
    if(eqa) wfcdf <- ecdf(wf[!is.na(wf)])
    if (cscale) {
        mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        on.exit(par(par.orig))
        w <- (3 + mar.orig[2]) * par("csi") * 2.54
        layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
        par(las = 1)
        mar <- mar.orig
        mar[4] <- mar[2]
        mar[2] <- 1
        par(mar = mar)
        levels <- seq(zlim[1], zlim[2], length=length(col)+1)
        plot.new()
        plot.window(xlim=c(0,1),ylim=range(levels), xaxs="i", yaxs="i")
        if (eqa) {
            vcol <- col[round((length(col)-1)*wfcdf(seq(zlim[1],zlim[2],length=length(col)))+1)]
        } else vcol <- col
	rect(0, levels[-length(levels)], 1, levels[-1], col=vcol, density=NA)
        axis(4)
        box()
        mar <- mar.orig
        mar[4] <- 0
        par(mar = mar)
    }
    if (is.null(cp)) {
        axis1 <- 1:nr
        axis2 <- 1:nc
    } else {
        axis1 <- ((1:nr)-cp$xc)/cp$rx
        axis2 <- ((1:nc)-cp$yc)/cp$ry
    }
    if (eqa) {
        iwf <- wfcdf(wf[!is.na(wf)])
        iwfm <- wf
        iwfm[!is.na(iwfm)] <- iwf
        zlim <- wfcdf(zlim)
        col1 <- round((length(col)-1)*zlim[1]+1)
        col2 <- round((length(col)-1)*zlim[2]+1)
        image(axis1, axis2, iwfm, zlim=zlim, asp=1, col=col[col1:col2], 
              xlab="X", ylab="Y", useRaster=TRUE, ...)
    } else {
        image(axis1, axis2, wf, zlim=zlim, asp=1, col=col, 
                  xlab="X", ylab="Y", useRaster=TRUE, ...)
    }
    if (addContours) contour(axis1, axis2, wf, add=TRUE)
}

## RGL animated 3D plot

## administrative stuff needed to make this a method for class "pupil"

wf3d <- function(wf, ...) UseMethod("wf3d", wf)

col3d <- function(wf, surf.col=topo.colors(256), zlim = NULL, eqa=FALSE) {
	if (is.null(zlim)) zlim <- range(wf, na.rm=TRUE)
    if (eqa) {
        wfcdf <- ecdf(wf[!is.na(wf)])
        iwf <- wfcdf(wf[!is.na(wf)])
        wf[!is.na(wf)] <- iwf
        zlim <- wfcdf(zlim)
    }
    surf.col[(length(surf.col)-1)*(wf-zlim[1])/(zlim[2]-zlim[1])+1]
}

wf3d.pupil <- function(wf, cp=NULL, zoom.wf=1, surf.col=topo.colors(256), bg.col="black",
                eqa=FALSE) {
    require(rgl)
    zlim <- range(wf, na.rm=TRUE)
    col <- col3d(wf, surf.col, zlim, eqa)
    if (is.null(cp)) {
        axis1 <- seq(-1, 1, length=nrow(wf))
        axis2 <- seq(-1, 1, length=ncol(wf))*(ncol(wf)/nrow(wf))
    } else {
        axis1 <- ((1:nrow(wf))-cp$xc)/cp$rx
        axis2 <- ((1:ncol(wf))-cp$yc)/cp$ry
    }

    rgl.bg(sphere=FALSE, fogtype="exp2", color=bg.col)
    rgl.surface(-axis1, axis2, wf*zoom.wf, color=col, shininess=100)
    rgl.lines(c(-1,-1.25)*max(axis1),c(0,0),c(0,0),color="red", lit=FALSE)
    rgl.lines(c(0,0),c(0,0),c(1,1.25)*max(axis2),color="red", lit=FALSE)
    rgl.texts(-1.3*max(axis1),0,0, "X", color="red")
    rgl.texts(0,0,1.3*max(axis2), "Y", color="red")
}

## polar coordinates inside the unit radius pupil

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
