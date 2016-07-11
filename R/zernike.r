## Assorted routines for manipulation of
## Zernike polynomials, interferogram analysis,
## and other forms of optical testing

## Author: M.L. Peck (mpeck1@ix.netcom.com)
## Language: R (http://www.r-project.org/)
## Copyright (c) 2004-2012, M.L. Peck
## Last mod: November 2012

## This is free software and is offered with no warranty whatsoever.

##utility function returns true if n-m is odd

.odd <- function(n,m) {
	if (((n-m)%%2) == 0) return(FALSE)
	else return(TRUE)
}

##radial zernike polynomials - iterative rewrite

rzernike <- function(rho, n, m) {
	if ((n<0) || (m<0) || (n<m) || (.odd(n,m))) return(0) #this is an error condition
	if (n==m) return(rho^n)
	j <- m
	rj <- rho^j
	r2 <- rho^2
	rjm2 <- 0
	for (j in seq(m,(n-2),by=2)) {
		c2 <- 4*(j+1)*r2-(j-m+2)^2/(j+2)
		c4 <- 0
		if (j != 0) {
			c2 <- c2 - (j+m)^2/j
			c4 <- (m^2 - j^2)/j
		}
		rjp2 <- (j+2)/((j+2)^2-m^2)*(c2 * rj + c4 * rjm2)
		rjm2 <- rj
		rj <- rjp2
	}
	return (rj)
}

##derivative (wrt r) of zernike polynomials
## Iterative version of recurrence relation
## given in Eqn. 13 of Mathworld Zernike article.

drzernike <- function(rho, n, m) {
	if ((n<0) || (m<0) || (n<m) || (.odd(n,m))) return(0)
	if ((n==0) && (m==0))return(0*rho)
	if (n==m) return(n*rho^(n-1))
	if (m==0) dr <- 0*rho else dr <- m*rho^(m-1)
	for (nn in seq(m+2, n, by=2))
		dr <- dr + nn*(rzernike(rho,nn-1,abs(m-1))+rzernike(rho,nn-1,m+1))
	return(dr)
}


### Zernike polynomial

Zernike <- function(rho, theta, n, m, t) {
	return( sqrt(n+1) * rzernike(rho, n, m) * switch(t, n=1, c=sqrt(2)*cos(m*theta),
	 s=sqrt(2)*sin(m*theta)))
}

### Radial derivative of zernike polynomial

DZernike <- function(rho, theta, n, m, t) {
	return( sqrt(n+1) * drzernike(rho, n, m) * switch(t, n=1, c=sqrt(2)*cos(m*theta),
	 s=sqrt(2)*sin(m*theta)))
}

### Tangential derivative of zernike polynomial

DTZernike <- function(rho, theta, n, m, t) {
	if (m == 0) return (rho*0) else
	return( sqrt(n+1) * m * rzernike(rho, n, m) * switch(t, n=0, c= -sqrt(2)*sin(m*theta), 
	 s = sqrt(2)*cos(m*theta)))
}



## create a unit aperture inside a matrix of arbitrary size
## note the parameter cp is a list with components (xc, yc, rx, ry, obstruct) as returned by pupil.pars()

pupil <- function(zcoef=NULL, zlist=makezlist(), phi=0, piston=0,
	nrow=256, ncol=nrow, cp=list(xc=128,yc=128,rx=127,ry=127,obstruct=0),
    obstruct=NULL) {

xs <- ((1:nrow)-cp$xc)/cp$rx
ys <- ((1:ncol)-cp$yc)/cp$ry

rhol <- function(x,y) {
 return(sqrt(x^2+y^2))
}

rho <- outer(xs,ys,rhol)
rho[rho>1] <- NA
if (is.null(obstruct))
	rho[rho<cp$obstruct] <- NA
else rho[rho<obstruct] <- NA

theta <- outer(xs, ys, function(x,y) atan2(y, x))
rho.v <- rho[!is.na(rho)]
theta.v <- theta[!is.na(rho)]
wf.v <- numeric(length(rho.v))
if (!is.null(zcoef)) wf.v <- zpm(rho.v, theta.v, phi, maxorder= max(zlist$n)) %*% c(0, zcoef)
wf <- matrix(NA, nrow=nrow, ncol=ncol)
wf[!is.na(rho)] <- wf.v
wf <- wf+piston
class(wf) <- "pupil"
return(wf)
}

pupil.arb <- function(zcoef=NULL, zlist=makezlist(), phi=0, piston=0,
	nrow=256, ncol=nrow, cp=list(xc=128,yc=128,rx=127,ry=127,obstruct=0),
    obstruct=NULL) {

xs <- ((1:nrow)-cp$xc)/cp$rx
ys <- ((1:ncol)-cp$yc)/cp$ry

rhol <- function(x,y) {
 return(sqrt(x^2+y^2))
}

rho <- outer(xs,ys,rhol)
rho[rho>1] <- NA
if (is.null(obstruct))
	rho[rho<cp$obstruct] <- NA
else rho[rho<obstruct] <- NA

theta <- outer(xs, ys, function(x,y) atan2(y, x))
rho.v <- rho[!is.na(rho)]
theta.v <- theta[!is.na(rho)]
wf.v <- numeric(length(rho.v))
if (!is.null(zcoef)) wf.v <- zpm.arb(rho.v, theta.v, phi, zlist=zlist) %*% zcoef
wf <- matrix(NA, nrow=nrow, ncol=ncol)
wf[!is.na(rho)] <- wf.v
wf <- wf+piston
class(wf) <- "pupil"
return(wf)
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
	return(exp(-(2*pi*rms)^2))
}

## A rainbow that may be better than R's rainbow() color palette

rygcb <- colorRampPalette(c("red", "yellow", "green", "cyan", "blue"), space = "Lab",
            interpolate="spline")
            
## Yet another rainbow made up of additive and subtractive primaries

rygcbm <- colorRampPalette(c("red", "yellow", "green", "cyan", "blue", "magenta"),
			space="Lab", interpolate="spline")

## a plot method for pupils

plot.pupil <- function(wf, cp=NULL, col=topo.colors(256), addContours=TRUE, cscale=FALSE, eqa=FALSE, zlim=NULL, ...) {
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
		image(axis1, axis2, iwfm, zlim=zlim, asp=1, col=col[col1:col2], xlab="X", ylab="Y", useRaster=TRUE, ...)
	} else image(axis1, axis2, wf, zlim=zlim, asp=1, col=col, xlab="X", ylab="Y", useRaster=TRUE, ...)
	if (addContours) contour(axis1, axis2, wf, add=TRUE)
}

## comparison plot of n wavefront estimates

plotn <- function(...,
  labels=NULL, addContours=FALSE,
  wftype="net", col=rygcb(400), qt=c(.01, .99)) {
	fits <- list(...)
	nwf <- length(fits)
	wfget <- paste("wf",wftype,sep=".")
	wfi <- get(wfget, pos=fits[[1]])
	nr <- nrow(wfi)
	nc <- ncol(wfi)
	wfs <- array(0, dim=c(nr,nc,nwf))
	wfs[,,1] <- wfi
	for (i in 2:nwf) wfs[,,i] <- get(wfget,pos=fits[[i]])
	rm(fits)
	if(is.null(labels)) labels=1:nwf

	rdiff <- matrix(0,nwf*(nwf-1)/2,2)
	k <- 1
	for (i in 1:(nwf-1)) {
		for (j in (i+1):nwf) {
			rdiff[k,] <- quantile(wfs[,,i]-wfs[,,j], probs=qt, na.rm=TRUE)
			k <- k+1
		}
	}
	zlim <- range(rdiff)
	if (tolower(.Platform$OS.type)=="windows") windows(5*nwf, 5*nwf) else X11(5*nwf,5*nwf)
	par(mar=c(0,0,2,0))
	split.screen(figs=c(nwf,nwf))
	for (i in 1:nwf) {
	  screen(i)
	  wfi <- wfs[,,i]
	  class(wfi) <- "pupil"
	  plot(wfi, col=col, addContours=addContours, eqa=(wftype=="net"), axes=FALSE)
	  title(main=paste(labels[i],"rms=",format(pupilrms(wfi),digits=3)))
	}
	scr <- nwf+1
	for (row in 1:(nwf-1)) {
	  screen(scr)
	  plot(0:1,0:1,type="n",axes=FALSE)
	  text(0.5,0.5,labels[row])
	  scr <- scr+row
	  for (i in (row+1):nwf) {
		screen(scr)
		image(1:nr,1:nc,wfs[,,row]-wfs[,,i],col=grey256,asp=1,axes=FALSE,zlim=zlim,useRaster=TRUE)
		title(main=paste("rms diff =",format(pupilrms(wfs[,,row]-wfs[,,i]), digits=3)))
		scr <- scr+1
	  }
	}
	close.screen(all.screens=T)
}

## plot cross sections through a wavefront map

plotxs <- function(wf, cp, theta0=0, ylim=NULL, N=4, n=101,
	col0="black", col="gray", lty=2) {
  theta0 <- theta0*pi/180
  ixy <- seq(-1,1,length=n)
  theta <- (0:(N-1))*pi/N
  ix <- cp$rx*ixy
  iy <- cp$ry*ixy
  if (is.null(ylim)) ylim <- range(wf, na.rm=TRUE)
  plot(range(ixy),ylim, type="n", xlab="rho", ylab="Height")
  for (i in 1:N) {
	iix <- round(cp$xc+cos(theta[i])*ix)
	iiy <- round(cp$yc+sin(theta[i])*iy)
	points(ixy, diag(wf[iix,iiy]), type="l", col=col, lty=lty)
  }
  iix <- round(cp$xc+cos(theta0)*ix)
  iiy <- round(cp$yc+sin(theta0)*iy)
  points(ixy, diag(wf[iix,iiy]), type="l", col=col0, lty=lty)
}



## a summary method for pupils

summary.pupil <- function(wf) {
	cat("Size:   ",nrow(wf), "x",ncol(wf),"\n")
	cat("RMS    =", format(pupilrms(wf), digits=3), "\n")
	cat("P-V    =", format(pupilpv(wf), digits=3), "\n")
	cat("Strehl =", format(strehlratio(pupilrms(wf)), digits=3), "\n")
}



## make a list of all orders up to maxorder

makezlist <- function(minorder=2, maxorder=14) {
	n <- numeric()
	m <- numeric()
	t <- character(length=0)

	for (order in seq(minorder, maxorder, by=2)) {
		mmax <- order/2
		mtemp <- numeric()
		for (j in mmax:1) mtemp <- c(mtemp, c(j, j))
		mtemp <- c(mtemp, 0)
		n <- c(n, order-mtemp)
		m <- c(m, mtemp)
		t <- c(t, rep(c("c", "s"), mmax), "n")
	}
return (list(n=n, m=m, t=t))
}

## Augmented Fringe set

zlist.fr <- makezlist(2,12)


## Vector of factors from conversion between "normalized" and conventional Zernike definitions

zmult <- function(zlist = makezlist()) {
	mult <- sqrt(zlist$n+1)
	mult[zlist$m > 0] <- sqrt(2)*mult[zlist$m > 0]
return(mult)
}



## create a matrix of Zernike polynomial values


zpm.arb <- function(rho, theta, phi=0, zlist=makezlist()) {
	zm <- matrix(0, nrow=length(rho), ncol=length(zlist$n))
	for (i in (1:length(zlist$n))) {
		zm[,i] <- Zernike(rho,theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
	}
	colnames(zm) <- colnames(zm, do.NULL=FALSE, prefix="Z")
	return(zm)
}

## a faster Zernike matrix fill

zpm <- function(rho, theta, phi=0, maxorder = 14) {
  nr <- length(rho)
  ncol <- (maxorder/2+1)^2
  zm <- numeric(nr*ncol)
  if (phi != 0) theta <- theta - pi*phi/180
  temp <- .C("zpmatrix", as.double(rho), as.double(theta), as.integer(nr),
	as.integer(maxorder), as.double(zm), package="zernike")
  zm <- temp[[5]]
  dim(zm) <- c(nr, ncol)
  colnames(zm) <- paste("Z",0:(ncol-1), sep="")
  return(zm)
}


## create a matrix of radial derivatives of Zernikes

filldzm <- function(rho, theta, phi=0, zlist=makezlist()) {
	dzm <- matrix(0, nrow=length(rho), ncol=length(zlist$n))
	for (i in (1:length(zlist$n))) {
		dzm[,i] <- DZernike(rho,theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
	}
	colnames(dzm) <- colnames(dzm, do.NULL=FALSE, prefix="DZ")
	return(dzm)
}

## create a matrix of gradient in polar coordinates of zernikes.
## radial derivative is on top, followed by tangential derivative.

fillgradientzm <- function(rho, theta, phi=0, zlist=makezlist()) {
	nr <- length(rho)
	drzm <- matrix(0, nrow=nr, ncol=length(zlist$n))
	dtzm <- matrix(0, nrow=nr, ncol=length(zlist$n))
	for (i in (1:length(zlist$n))) {
		drzm[,i] <- DZernike(rho, theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])
		dtzm[,i] <- DTZernike(rho, theta-pi*phi/180, zlist$n[i], zlist$m[i], zlist$t[i])/rho
		dtzm[rho==0,i] <- 0
	}
	gradzm <- rbind(drzm, dtzm)
	colnames(gradzm) <- colnames(gradzm, do.NULL=FALSE, prefix="Z")
	return(gradzm)
}

## fit zernikes to data

fitzernikes <- function(wf, rho, theta, phi=0, maxorder = 14 , uselm=FALSE) {
	zm <- zpm(rho, theta, phi=phi, maxorder=maxorder)
	zm.names <- colnames(zm)
	if (uselm) {
		fmla <- as.formula(paste("wf ~ -1 + ", paste(zm.names, collapse="+")))
		dataf <- data.frame(cbind(wf, zm))
		fit <- lm(fmla, data=dataf)
	}
	else {

		fit <- crossprod(solve(crossprod(zm)),crossprod(zm, wf))
	}
	return(fit)
}

## Zernike coefficients for a conic

zconic <- function(D, rc, b = -1, lambda = 1E-6, nmax = 6) {
	if ((nmax %%2) != 0) stop("nmax must be even")
	nterms <- nmax/2+1
	sa <- zc <- zs <- numeric(0)
	cz <- matrix(0, nterms, nterms)

	sd <- D/2
	na <- (sd/rc)^2
	sa[1] <- sa[2] <- 0
## 4th order conic term
	zs[3] <- (sd/rc)^3*sd/8/lambda
	zc[3] <- (1+b)*zs[3]
	sa[3] <- zc[3]-zs[3]
## binomial expansion of conic
	if (nterms >= 4) {
		for (k in 4:nterms) {
			mult <- (1-3/(2*(k-1)))*na
			zc[k] <- mult*(1+b)*zc[k-1]
			zs[k] <- mult*zs[k-1]
			sa[k] <- zc[k]-zs[k]
		}
	}	
## recurrence relation for coefficients of radial zernikes
	cz[1,] <- (-1)^(0:(nterms-1))
	cz[2,2] <- 2
	for (j in 3:nterms) {
		for (i in 2:j) {
			n <- 2*(j-1)
			cz[i, j] <- (4*(n-1)*cz[i-1, j-1]-2*(n-1)*cz[i,j-1]-(n-2)*cz[i,j-2])/n
		}
	}
	for (j in 1:nterms) {
		n <- 2*(j-1)
		cz[,j] <- cz[,j]*sqrt(n+1)
	}
	return(solve(cz,sa)[-(1:2)])
}

## returns cosine and sine components of Bath astigmatism

##
# D = diameter
# rc = radius of curvature
# s = beam separation
# lambda = source wavelength
# phi = angle from horizontal (degrees)
##

astig.bath <- function(D, rc, s, lambda=1e-6, phi=0) {
	astig.tot <- D^2*s^2/(32*sqrt(6)*lambda*rc^3)
	return(astig.tot*c(cos(pi*phi/90),sin(pi*phi/90)))
}



###########
## Various fun stuff
###########


## Star test simulator & support routines

grey256 <- gray(seq(0,1,length=256))
gray256 <- grey256

## calculate phase values for wavefront at wavelength lambda
## replaces NA values with 0

wftophase <- function(X, lambda = 1) {
	phi <- exp(2i*pi*X/lambda)
	phi[is.na(phi)] <- 0
	return(phi)
}

## puts matrix X into corner of npadded x npadded matrix
## padded with zeroes

padmatrix <- function(X, npad, fill=0) {
	nr <- nrow(X)
	nc <- ncol(X)
	xpad <- matrix(fill,npad,npad)
	xpad[1:nr,1:nc] <- X
	return(xpad)
}


## extract a matrix from the center of a larger matrix


submatrix <- function(X,size=255) {
	nr <- nrow(X)
	return(X[((nr-size)/2+1):((nr+size)/2),((nr-size)/2+1):((nr+size)/2)])
}

## shuffle quadrants of a 2d fft around to display as an image


fftshift <- function(X) {
	nr <- nrow(X)
	XS <- matrix(0,nr,nr)
	XS[1:(nr/2),1:(nr/2)] <- X[(nr/2+1):nr,(nr/2+1):nr]
	XS[(nr/2+1):nr,(nr/2+1):nr] <- X[1:(nr/2),1:(nr/2)]
	XS[(nr/2+1):nr,1:(nr/2)] <- X[1:(nr/2),(nr/2+1):nr]
	XS[1:(nr/2),(nr/2+1):nr] <- X[(nr/2+1):nr,1:(nr/2)]
	return(XS)
}

## computes & displays fraunhofer diffraction pattern
## & mtf for wavefront described in zernike coefficients zcoef

startest <- function(wf=NULL, zcoef=NULL, zlist=makezlist(), phi=0,
	lambda = 1, defocus=5,
	nrow = 255, ncol = nrow, cp = list(xc=128,yc=128,rx=127,ry=127,obstruct=0),
	obstruct=NULL, 
	npad = 4, 
	gamma=2, psfmag=2, displaymtf=TRUE, displaywf=FALSE) {

if (tolower(.Platform$OS.type) == "windows") {
  windows(width=15, height=5) } else {
  x11(width=15,height=5)
}
screens<- split.screen(c(1,3))

if (is.null(wf))
	wf <- pupil(zcoef=zcoef, zlist=zlist, phi=phi, nrow=nrow, ncol=ncol, cp=cp, obstruct=obstruct)
else {
	nrow <- nrow(wf)
	ncol <- ncol(wf)
}
wf.df <- pupil(zcoef=c(0,0,1), zlist=makezlist(2,2), nrow=nrow, ncol=ncol, cp=cp, obstruct=obstruct)
if (!is.null(obstruct)) wf[is.na(wf.df)] <- NA

lx <- round(2*cp$rx)+1
ly <- round(2*cp$ry)+1
npad <- npad * .up2(lx,ly)

phase <- wftophase(wf,lambda)
up <- Mod(fft(padmatrix(phase,npad)))
up <- up*up

otf <- fft(up, inverse=TRUE)/npad^2
otf <- otf[1:lx,1:ly]
mtf <- Re(otf)
mtf <- mtf/max(mtf)
freqx <- seq(0,1,length=lx)
freqy <- seq(0,1,length=ly)
mtfideal <- 2/pi*(acos(freqx)-freqx*sqrt(1-freqx^2))

nrpsf <- max(lx,ly)
psf <- submatrix(fftshift(up),floor(nrpsf/psfmag))
screen(screens[2])
image(psf^(1/gamma),col=grey256,asp=1,bty='n', axes=FALSE, useRaster=TRUE)
mtext("0")


if (defocus >5) nrpsf <- 2*nrpsf
if (defocus >15) nrpsf <- npad

phase <- wftophase(wf - defocus/3.46*lambda*wf.df, lambda)
up <- Mod(fft(padmatrix(phase,npad)))
up <- up*up
psf2 <- submatrix(fftshift(up),nrpsf)
screen(screens[1])
image(psf2^(1/gamma),col=grey256, asp=1,bty='n', axes=FALSE, useRaster=TRUE)
mtext(-defocus)

phase <- wftophase(wf + defocus/3.46*lambda*wf.df, lambda)
up <- Mod(fft(padmatrix(phase,npad)))
up <- up*up
psf2 <- submatrix(fftshift(up),nrpsf)
screen(screens[3])
image(psf2^(1/gamma),col=grey256,asp=1,bty='n',axes=FALSE, useRaster=TRUE)
mtext(defocus)
close.screen(all.screens=TRUE)


if (displaywf) {
	if (tolower(.Platform$OS.type) == "windows") windows() else x11()
	plot(wf/lambda)
}

if (displaymtf) {
	if (tolower(.Platform$OS.type) == "windows") windows() else x11()
	plot(freqx,mtf[1,],type="l",ylab="mtf",xlab="relative frequency")
	title(main='MTF vs. ideal')
	lines(freqy,mtf[,1])
	lines(freqx,mtfideal, lty=5)
	grid()
}
return(list(psf=psf, otf=otf, mtf=mtf))
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
	return(surf.col[(length(surf.col)-1)*(wf-zlim[1])/(zlim[2]-zlim[1])+1])
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

#	rgl.open()
	rgl.bg(sphere=FALSE, fogtype="exp2", color=bg.col)
	rgl.surface(-axis1, axis2, wf*zoom.wf, color=col, shininess=100)
	rgl.lines(c(-1,-1.25)*max(axis1),c(0,0),c(0,0),color="red", lit=FALSE)
	rgl.lines(c(0,0),c(0,0),c(1,1.25)*max(axis2),color="red", lit=FALSE)
	rgl.texts(-1.3*max(axis1),0,0, "X", color="red")
	rgl.texts(0,0,1.3*max(axis2), "Y", color="red")
}

## Kolmogorov turbulence

turbwf <- function(friedratio=1, zlist=makezlist(2,40), reps=1) {
	require(mvtnorm)
	dimz <- length(zlist$n)
	cova <- matrix(0, nrow=dimz, ncol=dimz)
	c0 <- friedratio^(5/3)*gamma(14/3)*(4.8*gamma(1.2))^(5/6)*(gamma(11/6))^2/(2^(8/3)*pi)
	for (i in 1:dimz) {
		for (j in i:dimz) {
			if ((zlist$m[i] == zlist$m[j]) && (zlist$t[i] == zlist$t[j])) {
				n <- zlist$n[i]
				np <- zlist$n[j]
				cova[i,j] <- (-1)^((n+np-2*zlist$m[i])/2) * sqrt((n+1)*(np+1)) *
					gamma((n+np-5/3)/2)/(gamma((np-n+17/3)/2)*gamma((n-np+17/3)/2)*
					gamma((n+np+23/3)/2))
				cova[j,i] <- cova[i,j]
			}
		}
	}
cova <- cova * c0/(4*pi^2)
zcoef.turb <- rmvnorm(reps, mean=rep(0,dimz), sigma=cova)
if (reps==1) zcoef.turb <- as.vector(zcoef.turb)
return(list(zcoef.turb=zcoef.turb, V=cova))
}


############## 
#Foucault test fits
#############

## modal fit to foucault data.


modalfit <- function(rho, f, rc, R, run=NULL, theta=rep(0, length(rho)),
    zlist=list(n=seq(4,12, by=2), m=rep(0, 5), t=rep("n", 5))) {
	if (is.null(run)) {
		nconst <- 1
		run <- as.factor(rep(1, length(rho)))
		mstring <- paste("DZ", 1:length(zlist$n), sep="")
	} else {
		nconst <- nlevels(run)
		mstring <- c('-1','run',paste("DZ", 1:length(zlist$n), sep=""))
	}
	eta <- f - (rho*R)^2/(2*rc)
	dzm <- -1e-6 * (rc/R)^2/rho * filldzm(rho, theta, phi=0, zlist=zlist)
	zm <- zpm.arb(rho, theta, phi=0, zlist=zlist)
	zernikes <- colnames(zm)
	dataframe <- data.frame(run, rho, f, eta, dzm, zm)
	fmla <- as.formula(paste("eta ~ ", paste(mstring, collapse= "+")))
	fit <- lm(fmla, data=dataframe)
	zcoef <- coef(fit)[-(1:nconst)]
	sigma.res <- summary(fit)$sigma
	cova <- sigma.res^2 * summary(fit)$cov.unscaled[-(1:nconst),-(1:nconst)]
	mse <- crossprod(zcoef)
	rms <- sqrt(mse)
	sigma.mse <- sqrt(4 * t(zcoef) %*% cova %*% zcoef + 2 * sum(diag(cova %*% cova)))
	sigma.rms <- .5 * sigma.mse/rms
	z <- data.matrix(dataframe[, zernikes])
	surface <- z %*% zcoef
	surface.se <- sqrt(diag( z %*% cova %*% t(z)))
	rho.fit <- seq(-1,1,length=201)
	zm <- zpm.arb(rho.fit, 0*rho.fit, zlist=zlist)
	surface.fit <- zm %*% zcoef
	return(list(rho.fit=rho.fit, surface.fit=surface.fit,
		lsfit=fit, zcoef=zcoef, rho.meas=rho, surface=surface, surface.se=surface.se, 
		mse=mse, sigma.mse=sigma.mse, rms=rms, sigma.rms=sigma.rms))
}



## Zonal fit to foucault data (integrate the moving source dif. eq.)

zonalfit <- function(rho, f, rc, R,
    zlist=list(n=seq(4,12, by=2), m=rep(0, 5), t=rep("n", 5))) {

	require(odesolve)

	couderpoints <- function(n, eps=0) {	### distribute points to approximate equal area over circle
		y <- sqrt(eps^2 + (1-eps^2)*(2*(1:n)-1)/(2*n))
		return(y)
	}


	ftrend.lo <- loess(f~I(rho^2), na.action=na.omit, control=loess.control(surface="direct"))
	foucault.difeq <- function(rho, z, parms) {
		rc <- parms[1]
		R <- parms[2]
		deriv <- R^2*rho/(rc+predict(ftrend.lo, newdata=data.frame(rho^2))-z)
		return(list(deriv, parms))
	}
	rho.fit <- couderpoints(501)
	z.fit <- lsoda(0, times=rho.fit, foucault.difeq, parms=c(rc, R), atol=1.e-8)[,2]

	#find the best fit paraboloid and subtract it off to get surface errors (scaled to nm)

	p.fit <- lm(z.fit ~ I(rho.fit^2))
	surface.fit <- 1e6*residuals(p.fit)

	## estimates of Zernike coefficients - note it only makes sense to use symmetrical ones
	zm <- zpm.arb(rho.fit, 0*rho.fit, zlist=zlist)
	zcoef <- coef(lm(surface.fit ~ zm))[-1]
	rms.surface <- sqrt(crossprod(surface.fit)/length(surface.fit))
	rho.fit <- c(-rev(rho.fit),rho.fit)
	surface.fit <- c(rev(surface.fit),surface.fit)
	return(list(rho.fit=rho.fit, surface.fit=surface.fit, zcoef=zcoef, rms.surf=rms.surface))
}

## Foucault test simulator

foucogram <- function(wf, edgex = 0, phradius = 0, slit=FALSE, pad=4, gamma=1, map =FALSE, lev=0.5) {

nr <- nrow(wf)
nc <- ncol(wf)
npad <- pad*nextn(max(nr,nc))
phi<-padmatrix(wftophase(wf),npad)
ca<-fftshift(fft(phi))
u <- npad/2 + 1 + edgex
ca[1:(u-phradius-1),]<-0
ca[u,] <- .5*ca[u,]
if (phradius>0) {
	if (slit) {
	for (i in 1:phradius) {
		f <- .5*(1+i/phradius)
		ca[(u+i),] <- f*ca[(u+i),]
		ca[(u-i),] <- (1-f)*ca[(u-i),]
	}
	} else
	for (i in 1:phradius) {
		f <- .5 + i*sqrt(phradius^2-i^2)/(pi*phradius^2) + asin(i/phradius)/pi
		ca[(u+i),] <- f*ca[(u+i),]
		ca[(u-i),] <- (1-f)*ca[(u-i),]
	}
}
ike <- Mod(fft(ca,inverse=T)[1:nr,1:nc])^2
image(ike^(1/gamma),col=grey256,asp=nc/nr,axes=FALSE, useRaster=TRUE)
if (map) {
	zmin <- floor(min(wf, na.rm=TRUE))
	zmax <- ceiling(max(wf, na.rm=TRUE))
	contour(wf,add=TRUE,levels=seq(zmin,zmax,by=lev),axes=FALSE,frame.plot=FALSE,col="red")
}
return(ike/max(ike))
}

## Synthetic interferogram

synth.interferogram <- function(wf=NULL, zcoef=NULL, zlist=NULL, 
  nr=nrow(wf), nc=ncol(wf), cp=NULL, phi=0, addzc=rep(0,4), fringescale=1, plots=TRUE) {
 	addzl <- makezlist(2,2)
	if (is.null(wf)) {
		if (is.null(cp)) {
			wf <- pupil(zcoef=zcoef, zlist=zlist, phi=phi)+
			  pupil(zcoef=addzc[-1], zlist=addzl, piston=addzc[1])
			nr <- nrow(wf)
			nc <- ncol(wf)
		}
		else {
			wf <- pupil(zcoef=zcoef, zlist=zlist, phi=phi, nrow=nr, ncol=nc, cp=cp) +
			        pupil(zcoef=addzc[-1], zlist=addzl, piston=addzc[1], nrow=nr, ncol=nc, cp=cp)
		}
	}
	else {
		if (is.null(cp))
			wf <- wf + pupil(zcoef=addzc[-1], zlist=addzl, piston=addzc[1], nrow=nr, ncol=nc)
		else
			wf <- wf + pupil(zcoef=addzc[-1], zlist=addzl, piston=addzc[1], nrow=nr, ncol=nc, cp=cp)
	}
	igram <- cos(2*pi*wf/fringescale)
	igram[is.na(igram)] <- 0
	class(igram) <- "pupil"
	if (plots) plot(igram, col=grey256, addContours=FALSE)
return(igram)
}

#########
## Rest are routines for FFT and PSI fringe analysis
#########

#########
## general utilities
#########

load.pgm <- function(files, imdiff=NULL) {
	require(pixmap)
	if (length(files)==1) {
		im1 <- t(attr(read.pnm(files), "grey"))
		im1 <- im1[,ncol(im1):1]
		return(im1)
	}
	if (!is.null(imdiff)) 
		im1 <- t(attr(read.pnm(files[imdiff]),"grey"))
	else
		im1 <- t(attr(read.pnm(files[1]), "grey"))
	nr <- nrow(im1)
	nc <- ncol(im1)
	im1 <- im1[,nc:1]
	if (!is.null(imdiff)) {
		images <- array(0, dim=c(nr, nc, length(files)-1))
		files <- files[-imdiff]
	}
	else
		images <- array(0, dim=c(nr, nc, length(files)))
	for (i in 1:length(files)) {
		images[,,i] <- t(attr(read.pnm(files[i]), "grey"))[,nc:1]
		if (!is.null(imdiff)) images[,,i] <- images[,,i]-im1
	}
	if (dim(images)[3]==1) images <- images[,,1]
	return(images)
}


## load image files (jpeg or tiff format) into an array

load.images <- function(files, names=files, channels=c(1,0,0), scale=1, FLIP=FALSE) {
	ext <- strsplit(names[1], ".", fixed=TRUE)[[1]]
	ext <- tolower(ext[length(ext)])
	fn <- switch(ext, jpg = , jpeg = readjpeg, tif = , tiff = readtiff)
	if (is.null(fn)) stop("Unrecognized file extension")
	readf <- function(fn, file, channels) {
	  im <- fn(file, channels)
	  if (scale != 1) im <- rescale(im, scale)
	  return(t(im[nrow(im):1,]))
	}
	if (length(files)==1) {
	   images <- readf(fn, files, channels)
	   ri <- range(images)
	   images <- (images-ri[1])/(ri[2]-ri[1])
	   if(FLIP) images <- images[nrow(images):1,]
	   return(images)
    }
	im1 <- readf(fn, files[1], channels)
	nr <- nrow(im1)
	nc <- ncol(im1)
	images <- array(0, dim=c(nr, nc, length(files)))
	images[,,1] <- im1
	for (i in 2:length(files)) {
		images[,,i] <- readf(fn, files[i], channels)
	}
	ri <- range(images)
	images <- (images-ri[1])/(ri[2]-ri[1])
	if(FLIP) images <- images[nr:1,,]
	return(images)
}

## rescale an image

rescale <- function(im, scale) {
  nr.new <- round(scale*nrow(im))
  nc.new <- round(scale*ncol(im))
  res <- .C("resize_image", as.double(im), 
			as.integer(ncol(im)), as.integer(nrow(im)),
			im.out = double(nr.new*nc.new),
			as.integer(nc.new), as.integer(nr.new),
			ret = integer(1), NAOK=TRUE, package="zernike")
  im.out = matrix(res$im.out, nr.new, nc.new)
  return(im.out)
}
  
## read a jpeg file

readjpeg <- function(filename, channels){
	res <- .C("read_jpg_img_info", as.character(filename),
	          width=integer(1), height=integer(1), depth=integer(1),
	          ret=integer(1), package="zernike")
	if (res$ret < 0)
	  stop(if (res$ret == -1) "Cannot open file." else "Internal error")
	width <- res$width
	height <- res$height
	depth <- res$depth
	res <- .C("read_jpg_img", as.character(filename),
			  as.integer(width), as.integer(height), as.integer(depth),
			  as.double(channels),
	          image=double(res$width * res$height),
	          ret=integer(1), package="zernike")
	image <- matrix(res$image, height, width)
	return(image)
}

## read a tiff file
readtiff <- function(filename, channels){
	res <- .C("read_tiff_img_info", as.character(filename),
	          width=integer(1), height=integer(1), depth=integer(1),
	          ret=integer(1), package="zernike")
	if (res$ret < 0)
		stop(if (res$ret == -1) "Cannot open file." else "Internal error")
	width <- res$width
	height <- res$height
	depth <- res$depth
	res <- .C("read_tiff_img", as.character(filename),
			  as.integer(width), as.integer(height), as.integer(depth),
			  as.double(channels),
	          image=double(res$width * res$height),
	          ret=integer(1), package="zernike")
	image <- matrix(res$image, height, width)
	return(image)
}

## crop an array

crop <- function(img, cp, npad=20) {
	nr <- dim(img)[1]
	nc <- dim(img)[2]
	xmin <- max(1, round(cp$xc-cp$rx-npad))
	xmax <- min(nr, round(cp$xc+cp$rx+npad))
	ymin <- max(1, round(cp$yc-cp$ry-npad))
	ymax <- min(nc, round(cp$yc+cp$rx+npad))
	cp.new <- cp
	cp.new$xc <- cp$xc-xmin+1
	cp.new$yc <- cp$yc-ymin+1
	if (length(dim(img)) > 2) img <- img[xmin:xmax,ymin:ymax,]
	  else img <- img[xmin:xmax,ymin:ymax]
	return(list(im=img, cp=cp.new))
}


pupil.pars <- function(im=NULL, obstructed=FALSE) {
	if (!is.null(im))
		image(1:nrow(im), 1:ncol(im), im, col=grey256, asp=1, xlab="", ylab="", useRaster=TRUE)
	cat("click on edge of pupil; right click to exit\n")
	flush.console()
	edge <- locator(type="p", col="green")
	x <- edge$x
	y <- edge$y
	el <- lm(x^2+y^2 ~ x + y + I(y^2))
	asp2 <- 1 - coef(el)[4]
	xc <- coef(el)[2]/2
	yc <- coef(el)[3]/(2*asp2)
	rx <- sqrt(coef(el)[1] + xc^2 + yc^2*asp2)
	ry <- rx / sqrt(asp2)
	if ((abs(rx-ry)/rx < 0.05) || (summary(el)$coefficients[4,4]>0.01)) {
		el <- update(el, ~ . - I(y^2))
		xc <- coef(el)[2]/2
		yc <- coef(el)[3]/2
		rx <- sqrt(coef(el)[1] + xc^2 + yc^2)
		ry <- rx
	}
	obstruct <- 0
	if (obstructed) {
		cat("click on edge of obstruction; right click to exit\n")
		flush.console()
		edge <- locator(type="p", col="green")
		is (!is.null(edge))
			obstruct <- sqrt(mean((edge$x-xc)^2+(edge$y-yc)^2*(rx/ry)^2))/rx
	}
	x.e <- seq(-rx,rx, length=101)
	y.e <- ry * sqrt(1 - (x.e/rx)^2)
	points(x.e+xc, y.e + yc, type='l', lty=1,col='green')
	y.e <- -y.e
	points(x.e+xc, y.e + yc, type='l', lty=1, col='green')
	if (obstruct > 0) {
		x.e <- obstruct * x.e
		y.e <- obstruct * y.e
		points(x.e+xc, y.e + yc, type='l', lty=1,col='green')
		y.e <- -y.e
		points(x.e+xc, y.e + yc, type='l', lty=1, col='green')
	}
	
	return(list(xc=xc, yc=yc, rx=rx, ry=ry, obstruct=obstruct))

}

## general purpose 2D convolution using FFT's.

## kern is the convolution kernel

convolve2d <- function(im, kern) {
	nr <- nrow(im)
	nc <- ncol(im)
	nrp <- nr + nrow(kern) - 1
	ncp <- nc + ncol(kern) - 1
	xs <- nrow(kern) %/% 2 + nrow(kern) %% 2
	ys <- ncol(kern) %/% 2 + ncol(kern) %% 2
	npad <- nextn(max(nrp, ncp))
	kern <- padmatrix(kern, npad=npad)
	im <- padmatrix(im, npad=npad)
	im.f <- Re(fft(fft(kern)*fft(im), inv=T))
	im.f <- im.f[xs:(nr+xs-1),ys:(nc+ys-1)]/(npad^2)
	return(im.f)
}


## Gaussian blur. fw is the standard deviation
## of the gaussian convolution kernel, in pixels.

gblur <- function(X, fw = 0, details=FALSE) {
  if (fw == 0) return(X)
  XP <- X
  XP[is.na(XP)] <- 0
  nr <- nrow(X)
  nc <- ncol(X)
  ksize <- max(ceiling(4 * fw), 3)
  if ((ksize %% 2) == 0) ksize <- ksize+1
  xc <- (ksize %/% 2) + 1
  xs <- ((1:ksize)-xc)/fw
  gkern <- outer(xs, xs, function(x,y) (x^2+y^2))
  gkern <- exp(-gkern/2)
  gkern <- round(gkern/min(gkern))
  npad <- nextn(max(nr,nc)+ksize-1)
  XP <- padmatrix(XP, npad)
  kernp <- padmatrix(gkern, npad)
  XP <- Re(fft(fft(XP)*fft(kernp), inv=TRUE))/(npad^2)/sum(gkern)
  XP <- XP[xc:(nr+xc-1), xc:(nc+xc-1)]
  XP[is.na(X)] <- NA
  if (details) 
	  return(list(gkern=gkern, X=XP))
  else return(XP)
}
 

 

## partial implementation of canny algorithm for edge detection
##
## arguments:
## im - the image to process
## fw - smoothing parameter for gaussian blur (pixels)
## qt - threshold for accepting a point as a candidate edge point.
## excl - number of pixels around edge to exclude. must be >= 1.
## plots - plot some intermediate results?
## details - return some extra data for testing purposes.

## there are many online references. Search "canny edge detection" or "canny algorithm".


circle.pars <- function(im, fw=2, qt=0.995, excl=5, refine=2,
  plots=TRUE, ask=TRUE, details=FALSE) {
	require(MASS)
	nr <- nrow(im)
	nc <- ncol(im)
	if (fw > 0) im <- gblur(im, fw)
	if (excl < 1) excl <- 1
	kern.sobel <- rbind(c(-1,-2,-1),c(0,0,0),c(1,2,1))
	gx <- convolve2d(im, kern.sobel)
	gy <- convolve2d(im, t(kern.sobel))
	modg <- sqrt(gx^2+gy^2)
	modg <- modg/max(modg)
	dirg <- atan2(gy, gx)
	
	# round off directions to 45 degrees
	
	dir <- dirg
	dir[abs(dirg) < 7*pi/8] <- (round(dirg[abs(dirg) < 7*pi/8]*4/pi))
	dir[abs(dirg) >= 7*pi/8] <- 0
	dir[dir < 0] <- 4+dir[dir<0]
	
	# find local maxima
	
	maxima <- NULL
	
	for (i in 0:3) {
		ptsi <- which(dir==i, arr.ind=TRUE)
		ed <- which(ptsi[,1] <= excl)
		ed <- c(ed, which(ptsi[,1] >= nr-excl+1))
		ed <- c(ed, which(ptsi[,2] <= excl))
		ed <- c(ed, which(ptsi[,2] >= nc-excl+1))
		if (length(ed) > 0)
			ptsi <- ptsi[-ed,]
		nb <- switch(i+1, c(1,0), c(1,1), c(0,1), c(1, -1))
		ptsp <- t(t(ptsi)+nb)
		ptsm <- t(t(ptsi)-nb)
		lm <- which((modg[ptsi] > modg[ptsp]) & (modg[ptsi] > modg[ptsm]))
		maxima <- rbind(maxima, ptsi[lm,])
	}
	
	# the thinned edges consist of the local maxima in the direction of gradient
	
	thin <- matrix(0, nr, nc)
	thin[maxima] <- modg[maxima]
	if (plots) image(1:nr, 1:nc, thin, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE)
		
	# the rest differs from published algorithm description. I'm just picking
	# edge candidates from top qt %-ile, and feeding them to 
	# lqs -- "robust" least squares routine. This is basically what I did
	# before.
	
	# future: Hough transform?? Probably not.
	
	ec <- which(modg[maxima] >= quantile(modg[maxima], probs=qt))
	maxima <- maxima[ec,]
	x <- maxima[,1]
	y <- maxima[,2]
	r2 <- x^2 + y^2
	ecm <- lqs(r2~x+y)
	xc <- coef(ecm)[2]/2
	yc <- coef(ecm)[3]/2
	rxy <- sqrt(coef(ecm)[1]+xc^2+yc^2)
	if (refine > 0) {
		xr <- (1:nr)-xc
		yr <- (1:nc)-yc
		rhod <- round(outer(xr, yr, function(x,y) sqrt(x^2+y^2)))
		edgec <- which(abs(rhod-rxy) <= refine, arr.ind=TRUE)
		ec <- which(thin[edgec] > 0)
		edgec <- edgec[ec,]
		x <- edgec[,1]
		y <- edgec[,2]
		ecm <- nls(0 ~ r^2-(x-xc)^2-(y-yc)^2, start=list(r=rxy, xc=xc, yc=yc),
		   weights=thin[edgec]^2)
		xc <- coef(ecm)[2]
		yc <- coef(ecm)[3]
		rxy <- coef(ecm)[1]
	}
		
	if (plots) {
		points(xc, yc, pch=20, col="red")
		if (ask) readline("Hit <enter> to continue ")
		points(x,y, pch=".", col="green")
		if (ask) readline("Hit <enter> to continue ")
		symbols(xc, yc, circles=rxy, inches=FALSE, add=TRUE, fg="red")
	}
	if (details)
		return(list(cp=list(xc=xc, yc=yc, rx=rxy, ry=rxy, obstruct=0),
		   modg=modg, thin=thin, maxima=maxima, finaledge=cbind(x,y), lsfit=ecm))
	else
		return(list(xc=xc, yc=yc, rx=rxy, ry=rxy, obstruct=0))
}

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
	return(list(rho=rho.mat, theta=theta.mat))
}


## Plot a complex matrix. 
## Maybe make cmat a class, so this becomes default plot routine?

plot.cmat <- function(X, fn="Mod", col=grey256, cp=NULL, zoom=1, gamma=1, ...) {
	  ## define some possibly useful functions. Note Mod2 produces power spectrum
	  
	logMod <- function(X) log(1+Mod(X))
	Mod2 <- function(X) Mod(X)^2
	logMod2 <- function(X) log(Mod2(X))
	nr <- nrow(X)
	nc <- ncol(X)
	xs <- seq(-nr/2,nr/2-1,length=nr)
	ys <- seq(-nc/2,nc/2-1,length=nc)
	if (!is.null(cp)) {
		xs <- (cp$rx/nr)*xs
		ys <- (cp$ry/nc)*ys
	}
	xsub <- max(1, floor(1+nr/2*(1-1/zoom))) : min(nr, ceiling(nr/2*(1+1/zoom)))
	ysub <- max(1, floor(1+nc/2*(1-1/zoom))) : min(nc, ceiling(nc/2*(1+1/zoom)))
	image(xs[xsub], ys[ysub], (eval(call(fn, X))[xsub,ysub])^(1/gamma), 
	  col=col, asp=1, xlab="k_x", ylab="k_y", useRaster=TRUE, ...)
}


pick.sidelobe <- function(imagedata, logm=FALSE, gamma=3) {
	imagedata <- imagedata - mean(imagedata)
	npad <- nextn(max(nrow(imagedata), ncol(imagedata)))
	im <- padmatrix(imagedata, npad=npad, fill=0)
	if (logm) fn <- "logMod2" else fn <- "Mod"
	im.fft <- fftshift(fft(im))
	plot.cmat(im.fft, fn=fn, gamma=gamma)
	cat("Click on the desired sidelobe peak\n")
	peak <- locator(n=1, type="p", col="red")
	cat("Click on background filter size\n")
	edge <- locator(n=1, type="n")
	hw <- round(sqrt((edge$x)^2+(edge$y)^2))
	symbols(0, 0, circles=hw, inches=FALSE, add=TRUE, fg="red")
	return(list(sl=c(round(peak$x), round(peak$y)), filter=hw))
}

## Difference that maintains original array size by appending a row

.fdiff <- function(X) rbind(diff(X), rep(NA,ncol(X)))

## Ghiglia and Pritt's wrap operator. Thanks to Steve K.

wrap <- function(phase) (phase+pi)%%(2*pi) - pi


## Calculate residue map
## A "positive" charged residue has value +1
## A "negative" charge has value -1

rmap <- function(phase, dx=NULL, dy=NULL, plot=FALSE, ...) {
  nr <- nrow(phase)
  nc <- ncol(phase)

  if (is.null(dx)) dx <- wrap(.fdiff(phase))
  if (is.null(dy)) dy <- wrap(t(.fdiff(t(phase))))
  d2y <- rbind(dy[1:(nr-1),]-dy[2:nr,],rep(NA,nc))
  d2x <- cbind(dx[,2:nc]-dx[,1:(nc-1)],rep(NA,nr))
  residues <- round((d2y+d2x)/(2*pi))
  if (plot) {
  	## positive charges
	chp <- which(residues ==1, arr.ind=TRUE)
	## negative charges
	chm <- which(residues == -1, arr.ind=TRUE)
	image(1:nr, 1:nc, phase, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE, ...)
	if(length(chp)>0) points(chp, col="green", pch=20)
	if(length(chm)>0) points(chm, col="red", pch=1)
	nr <- sum(abs(residues),na.rm=TRUE)
	return(nr)
  }
  else
  	return(residues)
}

## "Itoh's" method for phase unwrapping.

## integrate the wrapped phase differences
## This is called by brcutpuw, or it can be called
## directly by the user.

## Note: as of 15 May 2008 the unwrapper is written in C.
## It's fastest to bypass the computation of dx and dy in R
## and unwrap with a call to id_uw.
## netflowpuw (in lppuw) passes values of dx and dy with ucall==TRUE.
## brcutpuw passes dx and dy with ucall==FALSE.
## The call to id_dxy_uw is slower by about 80%.

idiffpuw <- function(phase, mask=phase, ucall=TRUE, dx=NULL, dy=NULL) {


	nr <- nrow(phase)
	nc <- ncol(phase)
	
	if (ucall) {
	   if (is.null(dx)) {
	       b <- .C("id_uw", as.integer(nr), as.integer(nc), as.double(phase/(2*pi)),
	           as.double(NA*phase), NAOK=TRUE, package="zernike")
	           puw <- matrix(b[[4]], nr, nc)
	   }
	   else {
	       b <- .C("id_dxy_uw", as.integer(nr), as.integer(nc), as.double(phase/(2*pi)),
	         as.double(mask), as.double(dx/(2*pi)), as.double(dy/(2*pi)),
	         as.double(NA*phase), integer(nr*nc), 
	         NAOK=TRUE, package="zernike")
	       puw <- matrix(b[[7]], nr, nc)
	   }
	   class(puw) <- "pupil"
	   return(puw)
	} else {
            if (is.null(dx)) dx <- wrap(.fdiff(phase))
            if (is.null(dy)) dy <- wrap(t(.fdiff(t(phase))))
	    b <- .C("id_dxy_uw", as.integer(nr), as.integer(nc), as.double(phase/(2*pi)), 
	      as.double(mask), as.double(dx/(2*pi)), as.double(dy/(2*pi)), 
	      as.double(NA*phase), integer(nr*nc), 
	      NAOK=TRUE, package="zernike")
	    puw <- matrix(b[[7]], nr, nc)
	    uw <- matrix(as.logical(b[[8]]), nr, nc)
	    class(puw) <- "pupil"
	    return(list(puw=2*pi*puw, uw=uw))
	}
}


## Quality guided unwrapper.
## This is mostly a wrapper for a call to the C function "q_uw".

qpuw <- function(phase, qual) {
	nr <- nrow(phase)
	nc <- ncol(phase)
	
	phase <- phase/(2*pi)
	qual[is.na(phase)] <- 0
	
	sol <- .C("q_uw", as.integer(nr), as.integer(nc), as.double(phase), as.double(qual),
		as.double(NA*phase), NAOK=TRUE, PACKAGE="zernike")
	puw <- matrix(sol[[5]], nr, nc)
	class(puw) <- "pupil"
	return(puw)
}


## All the high level routines calculate and display net, smoothed wavefronts extracted from
## raw phase estimates.

wf.net <- function(wf.raw, rho=NULL, theta=NULL, cp, sgs=1, zlist=makezlist(), zc0=c(1:3, 6:7), 
	satarget=c(0,0), astig.bath=c(0,0),
	uselm=FALSE, plots=TRUE) {
	nr <- nrow(wf.raw)
	nc <- ncol(wf.raw)
	subs <- matrix(FALSE,nr,nc)
	if (is.null(rho)) {
		prt <- pupil.rhotheta(nrow=nr, ncol=nc, cp=cp)
		rho <- prt$rho
		theta <- prt$theta
	}
	subs[seq(1,nr,by=sgs),seq(1,nc,by=sgs)] <- TRUE
	subs[is.na(rho)] <- FALSE
	wf.v <- wf.raw[subs]
	rho.v <- rho[subs]
	th.v <- theta[subs]

	## maybe a temporary fix. brcutpuw doesn't always unwrap completely
	
	rho.v <- rho.v[!is.na(wf.v)]
	th.v <- th.v[!is.na(wf.v)]
	wf.v <- wf.v[!is.na(wf.v)]
	fit <- fitzernikes(wf.v, rho.v, th.v, maxorder=max(zlist$n), uselm=uselm)
	if (uselm) cfit <- coef(fit) else cfit <- fit
	if (sign(cfit[9])*sign(satarget[1]) < 0) {
	  cfit <- -cfit
	  wf.raw <- -wf.raw
	  if (!uselm) fit <- -fit
	}
	zc.low <- rep(0,15)
	zc.low[zc0] <- cfit[zc0+1]
	zc.low[c(8,15)] <- zc.low[c(8,15)] + satarget
	zc.low[4:5] <- zc.low[4:5] + astig.bath
	wf.net <- wf.raw - pupil(zcoef=zc.low, zlist=makezlist(2,6), piston=cfit[1], 
	  nrow=nr, ncol=nc, cp=cp)
	if (plots) {
	  if (tolower(.Platform$OS.type) == "windows") {
	    windows(width=18, height=6) } else {
	    x11(width=18, height=6)
	  }
	  split.screen(figs=c(1,3))
	  screen(1)
	  plot(wf.net, cp=cp, addContours=FALSE)
	  mtext(paste("RMS = ", format(pupilrms(wf.net),digits=3)))
	}
	zcoef.net <- cfit[-1]
	zcoef.net[1:15] <- zcoef.net[1:15] - zc.low
	wf.smooth <- pupil(zcoef=zcoef.net, zlist=zlist, cp=cp, nrow=nr, ncol=nc)
	if (plots) {
	  screen(2)
	  plot(wf.smooth, cp=cp)
	  mtext(paste("RMS = ", format(pupilrms(wf.smooth),digits=3)))
	}
	wf.residual <- wf.net - wf.smooth
	if (plots) {
	  screen(3)
	  plot(wf.residual, cp=cp, col=grey256, addContours=FALSE)
	  mtext(paste("RMS = ", format(pupilrms(wf.residual),digits=3)))
	  close.screen(all.screens=TRUE)
	}
	return(list(wf.net=wf.net, wf.smooth=wf.smooth, wf.residual=wf.residual, 
	  fit=fit, zcoef.net=zcoef.net))
}

#############
## PSI routines
#############

## The following are low level routines to calculate the phase from various PSI algorithms.

## least squares fit of phase shifted interferograms with optional weights.

lspsi <- function(images, phases, wt=rep(1, length(phases))) {
	dims <- dim(images)
	X <- cbind(rep(1, length(phases)), cos(phases), sin(phases))
	X <- diag(wt) %*% X %*% solve(crossprod(diag(sqrt(wt)) %*% X))
	B <- array(matrix(images, dims[1]*dims[2], dims[3]) %*% X, dim=c(dims[1], dims[2], 3))
	phi <- atan2(-B[,,3], B[,,2])
	return(list(phi=phi, B=B))
}

## Vargas et al.'s (2011) Principal Components method

pcapsi <- function(im.mat, BGSUB=TRUE, diagpos=1) {
	if (BGSUB) im.mat <- im.mat-rowMeans(im.mat)

	# svd of the crossproduct is faster!
	
	svd.cp <- svd(crossprod(im.mat))
        if (diagpos == 1) {
            ph <- atan2(svd.cp$v[,2]*sqrt(svd.cp$d[2]),svd.cp$v[,1]*sqrt(svd.cp$d[1]))
            u <- im.mat %*% (svd.cp$u[,1:2] %*% diag(1/sqrt(svd.cp$d[1:2])))
        } else {
            ph <- atan2(svd.cp$v[,2],svd.cp$v[,1])
            u <- im.mat %*% svd.cp$u[,1:2]
        }
	ph <- wrap(ph-ph[1])
	phi <- atan2(-u[,2],u[,1])
	mod <- sqrt(u[,1]^2+u[,2]^2)
	r2 <- (svd.cp$d[1]+svd.cp$d[2])/sum(svd.cp$d)
	return(list(phi=phi, mod=mod/max(mod), phases=ph, snr=sqrt(r2/(1-r2)), eigen=svd.cp$d))
}

## The "advanced iterative algorithm" of Wang & Han (2004)

aiapsi <- function(im.mat, phases,
		maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE) {
	.meps <- sqrt(.Machine$double.eps)
	M <- nrow(im.mat)
	nf <- ncol(im.mat)
	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	sse <- numeric(maxiter+1)
	for (i in 1:maxiter) {

		# get the phase estimate from phase shifts
		
		phases.last <- phases
		SD <- svd(S)
		dp <- SD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
		sse[i] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
		if (i == 1) sse.1 <- sse[i]

		# toss the background and modulation portion of the ls solution
		
		phi <- atan2(-Phi[,3], Phi[,2])
		Phi <- cbind(rep(1,M), cos(phi), -sin(phi))

		# get phase shifts from the phase. Note use crossproduct for speed.
		
		PD <- svd(crossprod(Phi))
		dp <- PD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		S <- (PD$u %*% diag(dp) %*% t(PD$v)) %*% crossprod(Phi, im.mat)
		phases <- atan2(S[3,], S[2,])
		dphases <- sd(wrap(phases-phases.last))
		if (plotprogress) {
			if (i ==1)
				plot(1:maxiter, 1:maxiter, ylim=c(ptol, 1), type="n",
				  xlab="Iteration", ylab="", log="y")
			points(i, sse[i]/sse.1, pch=1)
			points(i, dphases, pch=2, col="green")
		}
		if ((trace > 0) && ((i-1)%%trace == 0)) {
		  cat(paste(i, ":", format(sse[i], digits=2), ":"), format(phases,digits=3), "\n")
		  flush.console()
		}
		if (dphases < ptol) break
		S <- rbind(rep(1,nf), cos(phases), sin(phases))
	}

	#final estimate of phase
	
	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	SD <- svd(S)
	dp <- SD$d
	dp[dp>.meps] <- 1/dp[dp>.meps]
	dp[dp <= .meps] <- 0
	Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
	sse[i+1] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
	phi <- atan2(-Phi[,3], Phi[,2])
	mod <- sqrt(Phi[,2]^2+Phi[,3]^2)
	return(list(phi=phi, mod=mod/max(mod), phases=phases, iter=i, sse=sse))
}


## A variation of Han & Kim's (1994) iterative algorithm

hkpsi <- function(im.mat, phases,
		maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE) {
	.meps <- sqrt(.Machine$double.eps)
	M <- nrow(im.mat)
	nf <- ncol(im.mat)
	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	sse <- numeric(maxiter+1)
	for (i in 1:maxiter) {
		phases.last <- phases
		
		# get the phase estimate from phase shifts
		
		SD <- svd(S)
		dp <- SD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
		sse[i] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
		if (i == 1) sse.1 <- sse[i]
		
		# get phase shifts from the phase. Note use crossproduct for speed.
		# Here is only difference from AIA - use full least squares estimate
		# at this step.
		
		PD <- svd(crossprod(Phi))
		dp <- PD$d
		dp[dp>.meps] <- 1/dp[dp>.meps]
		dp[dp <= .meps] <- 0
		S <- (PD$u %*% diag(dp) %*% t(PD$v)) %*% crossprod(Phi, im.mat)
		phases <- atan2(S[3,], S[2,])
		dphases <- sd(wrap(phases-phases.last))
		if (plotprogress) {
			if (i ==1)
				plot(1:maxiter, 1:maxiter, ylim=c(ptol, 1), type="n",
				  xlab="Iteration", ylab="", log="y")
			points(i, sse[i]/sse.1, pch=1)
			points(i, dphases, pch=2, col="green")
		}
		if ((trace > 0) && ((i-1)%%trace == 0)) {
		  cat(paste(i, ":", format(sse[i], digits=2), ":"), format(phases,digits=3), "\n")
		  flush.console()
		}
		if (dphases < ptol) break
		S <- rbind(rep(1,nf), cos(phases), sin(phases))
	}

	#final estimate of phase

	phases <- wrap(phases-phases[1])
	S <- rbind(rep(1,nf), cos(phases), sin(phases))
	SD <- svd(S)
	dp <- SD$d
	dp[dp>.meps] <- 1/dp[dp>.meps]
	dp[dp <= .meps] <- 0
	Phi <- tcrossprod(im.mat %*% (SD$v %*% diag(dp)), SD$u)
	sse[i+1] <- crossprod(as.vector(im.mat)-as.vector(Phi %*% S))
	phi <- atan2(-Phi[,3], Phi[,2])
	mod <- sqrt(Phi[,2]^2+Phi[,3]^2)
	return(list(phi=phi, mod=mod/max(mod), phases=phases, iter=i, sse=sse))
}

## PSI with variable tilt

tiltpsi <- function(im.mat, phases, x, y, tilts = NULL,
              nlpref=1, tlim=0.5, maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE) {
  if (nlpref>1) require(Rsolnp)
  M <- nrow(im.mat)
  nf <- ncol(im.mat)
  phases <- wrap(phases-phases[1])
  gotfirst <- FALSE

  # if no guess at tilts set them to 0 and calculate the phase the fast way
  
  if (is.null(tilts)) {
	tilts <- matrix(0, nf, 2)
	X <- cbind(rep(1,nf), cos(phases), sin(phases))
	X <- X %*% solve(crossprod(X))
	B <- im.mat %*% X
	gotfirst <- TRUE
  }
  tilts[,1] <- tilts[,1]-tilts[1,1]
  tilts[,2] <- tilts[,2]-tilts[1,2]

  # function that calculates phase given phase shifts and tilts
  
  phaseest <- function(im, phases, tilts) {
	M <- nrow(im)
	B <- numeric(3*M)
	lasol <- .C("pxls", as.integer(nf), as.integer(M), as.double(t(im)), as.double(phases),
	  as.double(tilts[,1]), as.double(tilts[,2]), as.double(x), as.double(y),
	  B=as.double(B), la.info=integer(1), package="zernike")
	B <- matrix(lasol$B, 3, M)
	if (lasol$la.info > 0) warning(paste("la.info =",lasol$la.info))
	return(t(B))
  }

  # this is the function minimized at the phase shift/tilt estimation step

  sseframe <- function(pt, im, phi, x, y, abar, bbar) {
	ph.xy <- pt[1] + 4*pi*(pt[2]*x+pt[3]*y)
	return(crossprod(im-abar-bbar*cos(phi+ph.xy)))
  }

  # gradient of the above function. used by nlminb

  gradsse <- function(pt, im, phi, x, y, abar, bbar) {
	ph.xy <- phi+pt[1]+4*pi*(pt[2]*x+pt[3]*y)
	cpt <- cos(ph.xy)
	spt <- sin(ph.xy)
	grads <- numeric(3)
	grads[1] <- 2*bbar*(-abar*sum(spt)+crossprod(im,spt)-bbar*crossprod(cpt,spt))
	grads[2] <- 8*pi*bbar*(-abar*crossprod(x,spt)+crossprod((im*x),spt)-bbar*crossprod((x*cpt),spt))
	grads[3] <- 8*pi*bbar*(-abar*crossprod(y,spt)+crossprod((im*y),spt)-bbar*crossprod((y*cpt),spt))
	return(grads)
  }


  pt.last <- cbind(phases, 2*pi*tilts)
  nltrace <- max(trace, 0)
  sse <- numeric(maxiter+1)
  sse.last <- 0

  for (i in 1:maxiter) {
	if ((i>1) || !gotfirst) B <- phaseest(im.mat, phases, tilts)

	# first column of B contains the background estimate; next two the components of phase
	# note: either mean or median seem to work here

	abar <- mean(B[,1])
	bbar <- mean(sqrt(B[,2]^2+B[,3]^2))
	phi <- atan2(-B[,3], B[,2])
	sse[i] <- 0

	# use one of two nonlinear optimizers for the tilt and phase shift
	# estimation step
	
	for (n in 1:nf) {
	  if (nlpref == 1) {
		smin <- nlminb(start=c(phases[n], tilts[n,]), objective=sseframe,
		  gradient=gradsse,
		  im=im.mat[,n], phi=phi, x=x, y=y, abar=abar, bbar=bbar,
		  control = list(trace=nltrace, abs.tol=1e-20),
		  lower=c(-2*pi, -tlim, -tlim), upper=c(+2*pi, tlim, tlim))
		sse[i] <- sse[i] + smin$objective
	  } else {
		smin <- solnp(pars=c(phases[n], tilts[n,]), fun=sseframe,
		  LB=c(-2*pi, -tlim, -tlim), UB=c(+2*pi, tlim, tlim),
		  control = list(trace=nltrace),
		  im=im.mat[,n], phi=phi, x=x, y=y, abar=abar, bbar=bbar)
		sse[i] <- sse[i] + smin$values[length(smin$values)]
	  }
	  phases[n] <- smin$par[1]
	  tilts[n,] <- smin$par[2:3]
	}
	if(smin$convergence > 0) warning("Convergence reported failed")

	# the tilts and phase shifts are offset from frame 1
	
	tilts[,1] <- tilts[,1]-tilts[1,1]
	tilts[,2] <- tilts[,2]-tilts[1,2]
	phases <- wrap(phases-phases[1])
	pt <- cbind(phases, 2*pi*tilts)
	dp <- sqrt(sum(diag(crossprod(pt-pt.last))))/(3*nf)
	pt.last <- pt
	if (plotprogress) {
		if (i ==1) {
			sse.1 <- sse[1]
			plot(1:maxiter, 1:maxiter, ylim=c(ptol, 1), type="n",
			  xlab="Iteration", ylab="", log="y")
		}
		points(i, sse[i]/sse.1, pch=1)
		points(i, dp, pch=2, col="green")
	}

	# this is slow, so it's best to print out some intermediate results so
	# the user knows something is happening.

	if (trace >= 0) print(paste("Iteration",i, "sse", format(sse[i], digits=2),
	  "delta sse",format((sse[i]-sse.last)/sse.last, digits=2),
	  "dp =", format(dp, digits=3)))
	if (dp < ptol) break
	sse.last <- sse[i]
  }
  B <- phaseest(im.mat, phases, tilts)
  phi <- atan2(-B[,3], B[,2])
  mod <- sqrt(B[,2]^2+B[,3]^2)
  return(list(phi=phi,mod=mod/max(mod),phases=phases,tilts=tilts,iter=i,sse=sse))

}


## High level routines

psifit <- function(images, phases, cp=NULL, wt=rep(1,length(phases)), 
		fringescale=1, zlist=makezlist(), zc0=c(1:3, 6:7), 
		satarget=c(0,0), astig.bath=c(0,0),
		puw.alg = "qual", uselm=FALSE, sgs=1, plots=TRUE, CROP=FALSE) {
	lsfit <- lspsi(images, phases, wt)
	phase.raw <- lsfit$phi
	mod.raw <- sqrt(lsfit$B[,,2]^2+lsfit$B[,,3]^2)
	mod.raw <- mod.raw/max(mod.raw)
	if (is.null(cp)) cp <- circle.pars(mod.raw, plot=plots, ask=FALSE)
	cp.orig <- cp
	if (CROP) {
	  mod.raw <- crop(mod.raw, cp)$im
	  phase.raw <- crop(phase.raw, cp)
	  cp <- phase.raw$cp
	  phase.raw <- phase.raw$im
	}
	nr <- nrow(phase.raw)
	nc <- ncol(phase.raw)
	prt <- pupil.rhotheta(nr, nc, cp)
	phase.raw[is.na(prt$rho)] <- NA
	class(phase.raw) <- "pupil"
	wf.raw <- switch(puw.alg,
	  qual = qpuw(phase.raw, mod.raw),
	  brcut = {require(lppuw); brcutpuw(phase.raw)},
	  lp = {require(lppuw); netflowpuw(phase.raw, mod.raw)})
	wf.raw <- fringescale*wf.raw
	class(wf.raw) <- "pupil"
	wf.nets <- wf.net(wf.raw, prt$rho, prt$theta, cp, sgs, zlist, zc0, satarget, astig.bath, uselm, plots)
	return(list(phase=phase.raw, mod=mod.raw, cp=cp,cp.orig=cp.orig,
	  wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth, 
	  wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net))
}


## PSI analysis by principal components (Vargas et al. 2011)

pcafit <- function(images, cp=NULL,
	BGSUB = TRUE, diagpos=1, REFINE = TRUE, 
        fringescale=1, zlist=makezlist(), zc0=c(1:3, 6:7), 
	satarget=c(0,0), astig.bath=c(0,0), 
	puw.alg="qual", uselm=FALSE, sgs=1, plots=TRUE, CROP=FALSE) {
	nr <- dim(images)[1]
	nc <- dim(images)[2]
	nf <- dim(images)[3]
	im.mat <- matrix(images,nr*nc,nf)
	if (!is.null(cp)) {
		prt <- pupil.rhotheta(nr,nc,cp)
		mask <- as.vector(prt$rho)
		im.mat <- im.mat[!is.na(mask),]
	}
	pcasol <- pcapsi(im.mat, BGSUB, diagpos)
	if (!is.null(cp)) {
	  phi <- matrix(NA, nr, nc)
	  mod <- matrix(0, nr, nc)
	  phi[!is.na(prt$rho)] <- pcasol$phi
	  mod[!is.na(prt$rho)] <- pcasol$mod
	} else {
	  phi <- matrix(pcasol$phi, nr, nc)
	  mod <- matrix(pcasol$mod, nr, nc)
	  cp <- circle.pars(mod, plot=plots, ask=FALSE)
	  if (REFINE) {
		prt <- pupil.rhotheta(nr,nc,cp)
		mask <- as.vector(prt$rho)
		im.mat <- im.mat[!is.na(mask),]
		pcasol <- pcapsi(im.mat, BGSUB, diagpos)
		phi <- matrix(NA, nr, nc)
		mod <- matrix(0, nr, nc)
		phi[!is.na(prt$rho)] <- pcasol$phi
		mod[!is.na(prt$rho)] <- pcasol$mod
	  }
	}
	class(phi) <- "pupil"
	cp.orig <- cp
	if (CROP) {
	  mod <- crop(mod, cp)$im
	  phi <- crop(phi, cp)
	  cp <- phi$cp
	  phi <- phi$im
	}
	nr <- nrow(phi)
	nc <- ncol(phi)
	prt <- pupil.rhotheta(nr,nc,cp)
	phi[is.na(prt$rho)] <- NA
    wf.raw <- switch(puw.alg,
      qual = qpuw(phi, mod),
      brcut = {require(lppuw); brcutpuw(phi)},
      lp = {require(lppuw); netflowpuw(phi, mod)})
	wf.raw <- fringescale*wf.raw
	class(wf.raw) <- "pupil"
	wf.nets <- wf.net(wf.raw, prt$rho, prt$theta, cp, sgs, zlist, zc0, satarget, 
	  astig.bath, uselm=uselm, plots=plots)
	return(list(phase=phi, mod=mod, phases=pcasol$phases, cp=cp, cp.orig=cp.orig,
	  wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth,
	  wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net, snr=pcasol$snr, eigen=pcasol$eigen))
}

## High level call for any of the iterative algorithms

itfit <- function(images, phases, cp=NULL,
		maxiter=20, ptol=0.001, trace=1, plotprogress=TRUE, REFINE=TRUE,
		tilts=NULL, tlim=0.5, nlpref=1,
		it.alg="aia",
		fringescale=1, zlist=makezlist(), zc0=c(1:3, 6:7),
		satarget=c(0,0), astig.bath=c(0,0),
		puw.alg = "qual", uselm=FALSE, sgs=1, plots=TRUE, CROP=FALSE) {
	nr <- dim(images)[1]
	nc <- dim(images)[2]
	nf <- dim(images)[3]
	im.mat <- matrix(images,nr*nc,nf)
	if (!is.null(cp)) {
		prt <- pupil.rhotheta(nr,nc,cp)
		mask <- as.vector(prt$rho)
		im.mat <- im.mat[!is.na(mask),]
	} else {
	  if (it.alg == "tilt") stop("Variable tilt algorithm must have non-null cp")
	}
	itsol <- switch(it.alg,
	  aia = aiapsi(im.mat, phases, maxiter=maxiter, ptol=ptol, trace=trace, plotprogress=plotprogress),
	  hk = hkpsi(im.mat, phases, maxiter=maxiter, ptol=ptol, trace=trace, plotprogress=plotprogress),
	  tilt = {x <- as.vector((prt$rho*cos(prt$theta))[!is.na(prt$rho)]);
			  y <- as.vector((prt$rho*sin(prt$theta))[!is.na(prt$rho)]);
			  tiltpsi(im.mat, phases, x, y, tilts=tilts, nlpref=nlpref, tlim=tlim,
			    maxiter=maxiter, ptol=ptol, trace=trace-1, plotprogress=plotprogress)
			 }
	)
	if (!is.null(cp)) {
	  phi <- matrix(NA, nr, nc)
	  mod <- matrix(0, nr, nc)
	  phi[!is.na(prt$rho)] <- itsol$phi
	  mod[!is.na(prt$rho)] <- itsol$mod
	} else {
	  phi <- matrix(itsol$phi, nr, nc)
	  mod <- matrix(itsol$mod, nr, nc)
	  cp <- circle.pars(mod, plot=plots, ask=FALSE)
	  if (REFINE) {
		prt <- pupil.rhotheta(nr,nc,cp)
		mask <- as.vector(prt$rho)
		im.mat <- im.mat[!is.na(mask),]
		itsol <- switch(it.alg,
	      aia = aiapsi(im.mat, phases, maxiter=maxiter, ptol=ptol, trace=trace, plotprogress=plotprogress),
	      hk = hkpsi(im.mat, phases, maxiter=maxiter, ptol=ptol, trace=trace, plotprogress=plotprogress)
		)
		phi <- matrix(NA, nr, nc)
		mod <- matrix(0, nr, nc)
		phi[!is.na(prt$rho)] <- itsol$phi
		mod[!is.na(prt$rho)] <- itsol$mod
	  }
	}
	cp.orig <- cp
	if (CROP) {
	  mod <- crop(mod, cp)$im
	  phi <- crop(phi, cp)
	  cp <- phi$cp
	  phi <- phi$im
	}
	nr <- nrow(phi)
	nc <- ncol(phi)
	prt <- pupil.rhotheta(nr,nc,cp)
	phi[is.na(prt$rho)] <- NA
    wf.raw <- switch(puw.alg,
      qual = qpuw(phi, mod),
      brcut = {require(lppuw); brcutpuw(phi)},
      lp = {require(lppuw); netflowpuw(phi, mod)})
	wf.raw <- fringescale*wf.raw
	class(wf.raw) <- "pupil"
	wf.nets <- wf.net(wf.raw, prt$rho, prt$theta, cp, sgs, zlist, zc0, satarget,
	  astig.bath, uselm=uselm, plots=plots)

	if(it.alg=="tilt") {
	  return(list(phase=phi, mod=mod, phases=itsol$phases, tilts=itsol$tilts, cp=cp, cp.orig=cp.orig,
		wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth,
		wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net,
		iter=itsol$iter, sse=itsol$sse))
	} else {
	  return(list(phase=phi, mod=mod, phases=itsol$phases, cp=cp, cp.orig=cp.orig,
		wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth,
		wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net,
		iter=itsol$iter, sse=itsol$sse))
	}
}


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



fftfit <- function(imagedata, cp=NULL, fringescale=1, sl=c(1,1), 
  filter=NULL, taper=2, 
  zlist=makezlist(), zc0=c(1:3, 6:7), 
  satarget=c(0,0), astig.bath=c(0,0),
  puw.alg = "qual", uselm=FALSE, sgs=2, plots=TRUE, CROP=FALSE) {
	
	nr <- nrow(imagedata)
	nc <- ncol(imagedata)
	npad <- nextn(max(nr,nc))
        if (!is.null(cp)) {
            prt <- pupil.rhotheta(nr, nc, cp)
            imagedata[is.na(prt$rho)] <- 0
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
	  if (sl[2] > 0) sfilter[, 1:(npad/2+1)] <- 0
	  else sfilter[, (npad/2+1):npad] <- 0
	} else
	if (sl[2] == 0) {
	  if (sl[1] > 0) sfilter[1:(npad/2+1),] <- 0
	  else sfilter[(npad/2+1):npad,] <- 0
	} 
	else 
	{
	  dydx <- -sl[1]/sl[2]
	  xcut <- (1:npad)-(npad/2+1)
	  ycut <- (npad/2+1) + round(dydx*xcut)
	  ycut[ycut<1] <- 1
	  ycut[ycut>npad] <- npad
	  if (sl[2] > 0)
		for (i in 1:npad) sfilter[i, 1:ycut[i]] <- 0
	  else
		for (i in 1:npad) sfilter[i, ycut[i]:npad] <- 0
	}
	im.fft <- im.fft*gblur(sfilter, fw=taper)
	sl.fft <- matrix(0, npad, npad)
	xmin <- max(1-sl[1], 1)
	xmax <- min(npad-sl[1], npad)
	ymin <- max(1-sl[2], 1)
	ymax <- min(npad-sl[2], npad)
	sl.fft[xmin:xmax, ymin:ymax] <- im.fft[(sl[1]+(xmin:xmax)),(sl[2]+(ymin:ymax))]
#	if (plots) plot.cmat(submatrix(sl.fft, size=npad/2))
	cphase <- fft(fftshift(sl.fft), inv=TRUE)[1:nr, 1:nc]
	phase.raw <- Arg(cphase)
	mod.raw <- Mod(cphase)
	mod.raw <- mod.raw/max(mod.raw)
	if (is.null(cp)) cp <- circle.pars(mod.raw, plot=plots, ask=FALSE)
	cp.orig <- cp
	if (CROP) {
	  mod.raw <- crop(mod.raw, cp)$im
	  phase.raw <- crop(phase.raw, cp)
	  cp <- phase.raw$cp
	  phase.raw <- phase.raw$im
	}
	nr <- nrow(phase.raw)
	nc <- ncol(phase.raw)
	prt <- pupil.rhotheta(nr, nc, cp)
	phase.raw[is.na(prt$rho)] <- NA
	class(phase.raw) <- "pupil"
    wf.raw <- switch(puw.alg,
      qual = qpuw(phase.raw, mod.raw),
      brcut = {require(lppuw); brcutpuw(phase.raw)},
      lp = {require(lppuw); netflowpuw(phase.raw, mod.raw)})
	wf.raw <- fringescale*wf.raw
	class(wf.raw) <- "pupil"
	wf.nets <- wf.net(wf.raw, prt$rho, prt$theta, cp, sgs, zlist, zc0, satarget, astig.bath, uselm, plots)
	return(list(phase=phase.raw, mod=mod.raw, cp=cp,cp.orig=cp.orig,
	  wf.net=wf.nets$wf.net, wf.smooth=wf.nets$wf.smooth, 
	  wf.residual=wf.nets$wf.residual, fit=wf.nets$fit, zcoef.net=wf.nets$zcoef.net))
}

## simple utility returns Euclidean length of a vector

hypot <- function(x) sqrt(crossprod(x))

## Zernike moments

zmoments <- function(zcoef, maxorder=14) {
    zcoef <- zcoef[-(1:3)]
    zlist <- makezlist(4, maxorder)
    nz <- length(zlist$n)
    sumstats <- NULL
    i <- 1
    repeat {
        if (i > nz) break
            if (zlist$m[i] == 0) {
                sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], zcoef[i]))
                i <- i+1
            }
            else {
                sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], hypot(zcoef[i:(i+1)])))
                i <- i+2
            }
    }
    sumstats <- data.frame(sumstats)
    names(sumstats) <- c("n", "m", "value")
    return(sumstats)
}

## these two functions provide rudimentary project management. The first adds
## Zernike coefficients and rotation angles from one or more fits to a matrix.
## The second separates polished in from instrumental aberrations (if possible)

addfit <- function(..., th=0, zcm=NULL, theta=numeric(0)) {
    fits <- list(...)
    nt <- length(th)
    theta <- c(theta, th*pi/180)
    nz <- length(fits[[1]]$zcoef.net)
    for (i in 1:nt) {
        zcm <- cbind(zcm, fits[[i]]$zcoef.net[4:nz])
    }
    return(list(zcm=zcm, theta=theta))
}

separate.wf <- function(zcm, theta, maxorder=14) {
    
    nt <- length(theta)
    zlist <- makezlist(4, maxorder)
    nz <- length(zlist$n)
    zcb <- matrix(0, nz, 4)
    sumstats <- NULL
    cx <- c(rep(1,nt),rep(0,nt))
    cy <- c(rep(0,nt),rep(1,nt))
    i <- 1
    repeat {
        if (i > nz) break
        if (zlist$m[i] == 0) {
            zcb[i,1] <- mean(zcm[i,])
            zcb[i,3] <- sd(zcm[i,])/sqrt(nt)
            sumstats <- rbind(sumstats, c(zlist$n[i], zlist$m[i], zcb[i,1],
                                rep(NA, 6)))
            i <- i+1
        } else {
            y <- c(zcm[i,],zcm[i+1,])
            rx <- c(cos(zlist$m[i]*theta),sin(zlist$m[i]*theta))
            ry <- c(-sin(zlist$m[i]*theta),cos(zlist$m[i]*theta))
            lsm <- lm(y ~ -1+rx+ry+cx+cy)
            cc <- coef(summary(lsm))[,1:2]
            zcb[i:(i+1),1] <- cc[1:2,1]
            zcb[i:(i+1),3] <- cc[1:2,2]
            if (nrow(cc)==4) {
                zcb[i:(i+1),2] <- cc[3:4,1]
                zcb[i:(i+1),4] <- cc[3:4,2]
            }
            lsm <- summary(lsm)
            sumstats <- rbind(sumstats,c(zlist$n[i],zlist$m[i],
                                hypot(zcb[i:(i+1),1]), hypot(zcb[i:(i+1),2]),
                                lsm$sigma, lsm$r.squared,lsm$fstatistic))
            i <- i+2
        }
    }
    zcb[is.na(zcb)] <- 0
    colnames(zcb) <- c("zc_mirror", "zc_inst", "se_zc_mirror", "se_zc_inst")
    colnames(sumstats)[1:7] <- c("n","m","rms_mirror", "rms_inst", "sigma","R2", "F")
    wf.mirror <- pupil.arb(zcoef=zcb[,1],zlist=zlist)
    wf.inst <- pupil.arb(zcoef=zcb[,2], zlist=zlist)
    return(list(zcb=data.frame(zcb),sumstats=data.frame(sumstats),
        wf.mirror=wf.mirror, wf.inst=wf.inst))
}
