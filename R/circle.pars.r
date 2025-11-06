## simple interactive edge finding using locator()

circle.man <- function(im=NULL, obstructed=FALSE) {
	if (!is.null(im))
		image(1:nrow(im), 1:ncol(im), im, col=grey256, asp=1, xlab="", ylab="", useRaster=TRUE)
	cat("click on edge of pupil; right click to exit\n")
	flush.console()
	edge <- locator(type="p", col="green")
	x <- edge$x
	y <- edge$y
	el <- lm(x^2+y^2 ~ x + y)
	xc <- coef(el)[2]/2
	yc <- coef(el)[3]/2
	rx <- sqrt(coef(el)[1] + xc^2 + yc^2)
	ry <- rx
	nlest <- nls(0 ~ r^2 - (x-xc)^2 - (y-yc)^2,
               start = list(r=rx, xc=xc, yc=yc))
  pars <- nlest$m$getPars()
  rx <- ry <- pars[1]
  xc <- pars[2]
  yc <- pars[3]
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
## sigma - smoothing parameter for gaussian blur (pixels)
## qt - threshold for accepting a point as a candidate edge point.
## excl - number of pixels around edge to exclude. must be >= 1.
## plots - plot some intermediate results?
## details - return some extra data for testing purposes.

## there are many online references. Search "canny edge detection" or "canny algorithm".


circle.auto <- function(img, sigma=3, qt=0.995, excl=5,
                        plots=TRUE, details=FALSE) {
  if (!require(robustbase)) {
    stop("Package robustbase is required for this function. Please install it from CRAN.")
  }

  nr <- nrow(img)
  nc <- ncol(img)

  dxy <- d_of_g(img, sigma)
  dirg <- atan2(dxy$Dy, dxy$Dx)
  modg <- dxy$Dxy

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
  # nlrob from package robustbase.


  ec <- which(modg[maxima] >= quantile(modg[maxima], probs=qt))
  maxima <- maxima[ec,]
  x <- maxima[,1]
  y <- maxima[,2]
  df <- data.frame(x=x, y=y)
  ecm <- nlrob(0 ~ r^2 - (x-xc)^2 - (y-yc)^2, data=df,
               start=list(r=min(nr,nc)/2, xc=nr/2, yc=nc/2))
  xc <- coef(ecm)['xc']
  yc <- coef(ecm)['yc']
  rxy <- coef(ecm)['r']

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

