##' Estimate parameters of a circle using Hough Circle Transform
##' 
##' Uses a portion of the canny algorithm to find candidate edge points
##' and the direction of the gradient at those points, then uses
##' Hough Circle Transform to estimate circle parameters.
##' 
##' @param im The image to find a circle in (a modulation estimate is best)
##' @param fw Size of Gaussian blur to smooth image
##' @param qt Threshold to accept strong edge candidate
##' @param excl Number of pixels to exclude around edge of frame as candidates
##' @param rmin Minimum circle radius
##' @param rmax Maximum circle radius
##' @param rstep step size in constructing lookup table
##' @param dtheta_max maximum assumed error in gradient direction
##' @param dtheta_step increment for dtheta
##' @param number of nearest neighbors for alternate calculation
##' @param plots plot?
##' @param details Return extra details?
##' 
##' @details
##' The Hough transform section first creates a lookup table of candidate
##' radii and center points, then for each candidate edge point calculates
##' potential centers along a fan of rays near the gradient direction.
##' An inner join then finds matches in the lookup table and increments
##' an accumulator vector. Highest vote at the end wins.
##' 
##' @return If details is FALSE a named list with the circle parameters
##' 
##' @note
##' This is experimental and can be very slow. A good guess for the radius is very helpful.
##' Experimental feature: find the nn nearest neighbors of the selected trio of parameters
##' and calculate a vote weighted mean. This is returned as \code{rxy_alt} if details is TRUE.
##' 
##' @seealso [circle.pars()], [pupil.pars()]
##' 
##' @examples
##' example("psifit", package="zernike", ask=FALSE)
##' X11()
##' cp2 <- circle.hough(tfit$mod, rmin=round(tfit$cp$rx)-10, rmax=round(tfit$cp$rx)+10)
circle.hough <- function(im, fw=2, qt=0.995, excl=5,
                         rmin=min(dim(im))/4, rmax=min(dim(im))/2, rstep=1,
                         dtheta_max = 0.5, dtheta_step=0.05,
                         nn=7, plots=TRUE, details=FALSE) {
  require(data.table)
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
  ec <- which(modg[maxima] >= quantile(modg[maxima], probs=qt))
  maxima <- maxima[ec,]
  dir_edges <- dirg[maxima]

  if (plots) {
    image(1:nr, 1:nc, thin, col=grey256, asp=1, xlab="X", ylab="Y", useRaster=TRUE)
    points(maxima[,1], maxima[,2], pch=20, col='green')
  }
  
  # Hough circle transform
  
  xc_min <- rmin
  xc_max <- nr-rmin
  yc_min <- rmin
  yc_max <- nc-rmin
  
  r_cand <- seq(rmin, rmax, by=rstep)
  xc_cand <- seq(xc_min, xc_max, by=rstep)
  yc_cand <- seq(yc_min, yc_max, by=rstep)
  rxy <- data.table(expand.grid(r_cand, xc_cand, yc_cand))
  rxy <- cbind(index=1:nrow(rxy), rxy)
  names(rxy) <- c("index", "r", "xc", "yc")
  dtheta <- seq(-dtheta_max, dtheta_max, by=dtheta_step)
  lt <- length(dtheta)
  acc <- numeric(nrow(rxy))
  
  for (i in 1:nrow(maxima)) {
    xc <- maxima[i,1] - as.vector(outer(r_cand, dir_edges[i]+dtheta, function(X,Y) X * cos(Y)))
    xc <- rstep * round(xc/rstep)
    yc <- maxima[i,2] - as.vector(outer(r_cand, dir_edges[i]+dtheta, function(X,Y) X * sin(Y)))
    yc <- rstep * round(yc/rstep)
    circles <- data.table(r=rep(r_cand, lt), xc=xc, yc=yc)
    matches <- dplyr::inner_join(rxy, circles, by=c("r", "xc", "yc"))
    acc[matches$index] = acc[matches$index] + 1
  }
  rxy_best <- rxy[which.max(acc),]
  cp <- list(xc = rxy_best$xc, yc=rxy_best$yc, rx=rxy_best$r, ry=rxy_best$r, obstruct=0)
  
  wtmean <- function(x, wt) sum(wt*x)/sum(wt)
  
  gets_vote <- as.vector(nabor::knn(rxy[,-1],rxy_best[,-1], k=nn)$nn.idx)
  rxy_alt <- list(r=wtmean(rxy$r[gets_vote],acc[gets_vote]), xc=wtmean(rxy$xc[gets_vote],acc[gets_vote]), 
                  yc=wtmean(rxy$yc[gets_vote],acc[gets_vote]))
  
  if (plots) {
    points(rxy_best$xc, rxy_best$yc, col="red", pch=20)
    symbols(rxy_best$xc, rxy_best$yc, circles=rxy_best$r, inches=FALSE, add=TRUE, fg="red")
  }
  if (details) {
    list(cp=cp, thin=thin, maxima=maxima, dir_edges=dir_edges, 
         rxy=rxy, acc=acc, centers=centers, rxy_best=rxy_best, rxy_alt=rxy_alt)
  } else {
    cp
  }
}
                        
                        
