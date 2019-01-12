# demo of PSI fitting process
# PSI data courtesy of Vladimir Galogaza.

# list of image files
# these are stored in zernike/psidata in the R library, so must expand path name

require(zernike)
fpath <- file.path(find.package(package="zernike"), "psidata")
files <- scan(file.path(fpath, "files.txt"), what="character")
for (i in 1:length(files)) files[i] <- file.path(fpath, files[i])

# load the images into an array

images <- load.images(files)

# parameters for this run

source(file.path(fpath, "parameters.txt"))

# phase shifts

phases <- wrap((0:(dim(images)[3]-1))/frames.per.cycle*2*pi)
phases <- switch(ps.dir, ccw = -phases, cw = phases, phases)

# target SA coefficients for numerical null.

sa.t <- sconic(diam,roc,lambda=wavelength)
zopt <- get_zoptions()
zopt$plots <- FALSE
zopt$satarget <- sa.t

# print it

sa.t #Note I correct 4th and 6th order Zernike SA terms

# display an interferogram

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
image(1:dim(images)[1], 1:dim(images)[2], images[,,1], col=grey256, asp=1,
 xlab="X", ylab="Y", useRaster=TRUE)
mtext("Sample Interferogram")

# this does the psi fit, and displays the results
cat("Performing PSI analysis. Please wait...\n")
flush.console()
psfit <- psifit(images, phases, psialg="ls", options=zopt)

# what's stored in psfit

attributes(psfit)

# Net zernike coefficients

psfit$zcoef.net

# Zernike coefficients in traditional scaling

psfit$fit[-1]*zmult(makezlist(2,14))

# phase map

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
mtext(paste("Wrapped phase map -", rmap(psfit$phi, plot=TRUE), "corrupted pixels", sep=" "))

# display the smoothed wavefront, and show some summary information

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(psfit$wf.smooth)
mtext(paste("Smoothed wavefront - RMS =",format(pupilrms(psfit$wf.smooth),digits=3)))
cat("Summary quality measures\n")
summary(psfit$wf.smooth)

# Summary is too optimistic since the test wavelength is in the red. Restate for 550nm.

cat("Summary quality measures at 550 nm\n")
summary(wavelength/550*psfit$wf.smooth)

# same for net, unsmoothed wavefront

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(psfit$wf.net,addContours=FALSE)
mtext(paste("Net wavefront - RMS =",format(pupilrms(psfit$wf.net),digits=3)))
cat("Summary quality measures for unfiltered wavefront\n")
summary(psfit$wf.net)

# and now residual map

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(psfit$wf.residual,addContours=FALSE,col=grey256)
mtext(paste("Residuals from fit - RMS =",format(pupilrms(psfit$wf.residual),digits=3)))

# do a 3D plot

# use dynamic version if rgl is available

wf.smooth <- psfit$wf.smooth
x0 <- round(psfit$cp$xc)-dim(images)[2]/2
wf.smooth <- wf.smooth[x0:(x0+dim(images)[2]-1),]
class(wf.smooth) <- "pupil"
if (any(.packages(all.available=TRUE) == "rgl")) wf3d(wf.smooth)

  
  # define a new function
  
perspwf <- function(wf, theta=0, phi=30, ...) {
	wfi <- (wf[-1, -1] + wf[-1, -ncol(wf)] + wf[-nrow(wf), -1] + wf[-nrow(wf), -ncol(wf)])/4
	zlim <- range(wfi, finite=TRUE)
	colors <- matrix("white", nrow=nrow(wf)-1, ncol=ncol(wf)-1)
	colors <- topo.colors(256)[255*(wfi-zlim[1])/(zlim[2]-zlim[1])+1]
	xyaxis <- seq(-1,1,length=nrow(wf))
	persp(xyaxis, xyaxis, wf, theta=theta, phi=phi, scale=FALSE,
	    col=colors, border=NA, shade=0.5,
	    box=FALSE, axes=FALSE, ...)
}
#persp plot
if (tolower(.Platform$OS.type) == "windows") windows() else x11()
perspwf(wf.smooth)
mtext("Perspective plot")

# calculate the expected Bath astigmatism and subtract it from wavefront

astig.bath <- diam^2*beam.sep^2/roc^3/(32*sqrt(6)*wavelength*1.e-6)
zcoef.net <- psfit$zcoef.net
zcoef.net[4] <- zcoef.net[4]-astig.bath
wf.smooth <- pupil(zcoef=zcoef.net, zlist=makezlist(2,14))

#plot it with a different color palette

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(wf.smooth, col=rygcb(256), cscale=TRUE)
mtext(paste("Bath astigmatism removed - RMS =",format(sqrt(crossprod(zcoef.net)),digits=3)))


# clean up

rm(x0, wf.smooth, zcoef.net)
