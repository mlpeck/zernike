## simple interactive edge finding using locator()

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
	
	list(xc=xc, yc=yc, rx=rx, ry=ry, obstruct=obstruct)

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


circle.pars <- function(im, fw=2, qt=0.995, excl=5,
                        plots=TRUE, details=FALSE) {
  nr <- nrow(im)
  nc <- ncol(im)
  if (fw > 0) {
    im <- gblur(im, fw)
  }
  if (excl < 1) {
    excl <- 1
  }
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
    if (length(ed) > 0) {
      ptsi <- ptsi[-ed,]
    }
    nb <- switch(i+1, c(1,0), c(1,1), c(0,1), c(1, -1))
    ptsp <- t(t(ptsi)+nb)
    ptsm <- t(t(ptsi)-nb)
    lm <- which((modg[ptsi] > modg[ptsp]) & (modg[ptsi] > modg[ptsm]))
    maxima <- rbind(maxima, ptsi[lm,])
  }
  
  # the thinned edges consist of the local maxima in the direction of gradient
  
  thin <- matrix(0, nr, nc)
  thin[maxima] <- modg[maxima]
  if (plots) {
    image(1:nr, 1:nc, thin, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE)
  }
  
  # the rest differs from published algorithm description. I'm just picking
  # edge candidates from top qt %-ile, and feeding them to 
  # nlrob from package robustbase if available or
  # lqs -- "robust" least squares routine. This is basically what I did
  # before.
  
  
  ec <- which(modg[maxima] >= quantile(modg[maxima], probs=qt))
  maxima <- maxima[ec,]
  x <- maxima[,1]
  y <- maxima[,2]
  if (require(robustbase)) {
    df <- data.frame(x=x, y=y)
    ecm <- nlrob(0 ~ r^2 - (x-xc)^2 - (y-yc)^2, data=df, 
                 start=list(r=min(nr,nc)/2, xc=nr/2, yc=nc/2))
    xc <- coef(ecm)['xc']
    yc <- coef(ecm)['yc']
    rxy <- coef(ecm)['r']
  } else {
    require(MASS)
    r2 <- x^2 + y^2
    ecm <- lqs(r2~x+y)
    xc <- coef(ecm)[2]/2
    yc <- coef(ecm)[3]/2
    rxy <- sqrt(coef(ecm)[1]+xc^2+yc^2)
  }
#   if (refine > 0) {
#     xr <- (1:nr)-xc
#     yr <- (1:nc)-yc
#     rhod <- round(outer(xr, yr, function(x,y) sqrt(x^2+y^2)))
#     edgec <- which(abs(rhod-rxy) <= refine, arr.ind=TRUE)
#     ec <- which(thin[edgec] > 0)
#     edgec <- edgec[ec,]
#     x <- edgec[,1]
#     y <- edgec[,2]
#     ecm <- nls(0 ~ r^2-(x-xc)^2-(y-yc)^2, start=list(r=rxy, xc=xc, yc=yc),
#                weights=thin[edgec]^2)
#     xc <- coef(ecm)['xc']
#     yc <- coef(ecm)['yc']
#     rxy <- coef(ecm)['r']
#   }
  
  if (plots) {
    points(xc, yc, pch=20, col="red")
    points(x,y, pch=20, col="green")
    symbols(xc, yc, circles=rxy, inches=FALSE, add=TRUE, fg="red")
  }
  if (details) {
    finaledge <- cbind(x,y)
    list(cp=list(xc=xc, yc=yc, rx=rxy, ry=rxy, obstruct=0),
         gx=gx, gy=gy, modg=modg, dirg=dirg,
         thin=thin, maxima=maxima, finaledge=finaledge, dir_edge=dirg[finaledge],
         lsfit=ecm)
  } else {
    list(xc=xc, yc=yc, rx=rxy, ry=rxy, obstruct=0)
  }
}
                        
                        
