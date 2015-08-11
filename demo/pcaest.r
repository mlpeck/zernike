# demo of PCA fitting
# phase shifted data courtesy of Vladimir Galogaza.

# list of image files
# these are stored in zernike/psidata in the R library, so must expand path name

require(zernike)

plotn <- function(wfs, legends=rep("",dim(wfs)[3]), wftype="net", col=rygcb(400), qt=c(.1, .9)) {
  nwf <- dim(wfs)[3]
  ncomp <- nwf*(nwf-1)/2
  k <- 1
  rdiff <- matrix(0, ncomp, 2)
  for (i in 1:(nwf-1)) {
	for (j in (i+1):nwf) {
	  rdiff[k,] <- quantile(wfs[,,i]-wfs[,,j], probs=qt, na.rm=TRUE)
	  k <- k+1
	}
  }
  zlim <- range(rdiff)
  if (tolower(.Platform$OS.type) == "windows") {
	windows(width=4*nwf,height=4*nwf)} else {
	x11(width=4*nwf, height=4*nwf)
  }
  par(mar=c(0,0,1,0))
  split.screen(figs=c(nwf,nwf))
  for (k in 1:nwf) {
	screen(k)
	plot.pupil(wfs[,,k], col=col, addContours=FALSE, eqa=(wftype=="net"), axes=F)
	mtext(paste(legends[k],"rms=",format(pupilrms(wfs[,,k]),digits=3)))
  }
  k <- k+1
  for (i in 1:(nwf-1)) {
	screen(k)
	plot(0:1,0:1,type="n",axes=F)
	text(0.5,0.5,legends[i])
	k <- k+i
	for (j in (i+1):nwf) {
	  screen(k)
	  plot.pupil(wfs[,,i]-wfs[,,j], col=grey256, addContours=FALSE, axes=FALSE, zlim=zlim)
	  mtext(paste("rms diff =",format(pupilrms(wfs[,,i]-wfs[,,j]), digits=3)))
	  k <- k + 1
	}
  }
  close.screen(all.screens=TRUE)
}

	  
fpath <- file.path(find.package(package="zernike"), "psidata")
files <- scan(file.path(fpath, "files_pca.txt"), what="character")
for (i in 1:length(files)) files[i] <- file.path(fpath, files[i])

# load the images into an array

images <- load.images(files)
nr <- dim(images)[1]
nc <- dim(images)[2]
nf <- dim(images)[3]

# parameters for this run

source(file.path(fpath, "parameters.txt"))

# phase shifts

phases <- wrap((0:(nf-1))/frames.per.cycle*2*pi)
phases <- switch(ps.dir, ccw = -phases, cw = phases, phases)

# target SA coefficients for numerical null.

sa.t <- 2*zconic(diam,roc,lambda=wavelength*1e-6)


cat("Performing PSI analysis. Please wait...\n")
flush.console()
system.time(psfit <- psifit(images,phases, satarget=sa.t,plots=FALSE))

cat("Zernike fit to wavefront from PSI analysis\n")
flush.console()

plot(psfit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(psfit$wf.smooth),digits=3)))

cat("Performing PCA analysis. Please wait...\n")
system.time(pcfit <- pcafit(images, cp=psfit$cp, satarget=sa.t, plots=FALSE))

cat("Wavefront from PCA routine\n")
flush.console()

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(pcfit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(pcfit$wf.smooth),digits=3)))

cat("Performing 'Advanced iterative algorithm'. Please wait...\n")
system.time(aifit <- itfit(images, phases=pcfit$phases, cp=psfit$cp, satarget=sa.t,
  plots=FALSE, trace=1))

cat("Wavefront from 'AIA'\n")
flush.console()

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(aifit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(aifit$wf.smooth),digits=3)))


# full intercomparison

wfnets <- array(0, dim=c(nr, nc, 3))
wfsmooth <- array(0, dim=c(nr,nc,3))

legends=c("baseline PSI", "PCA", "AIA")

wfnets[,,1] <- psfit$wf.net
wfnets[,,2] <- pcfit$wf.net
wfnets[,,3] <- aifit$wf.net

wfsmooth[,,1] <- psfit$wf.smooth
wfsmooth[,,2] <- pcfit$wf.smooth
wfsmooth[,,3] <- aifit$wf.smooth

plotn(wfnets,legends=legends)
plotn(wfsmooth,legends=legends,wftype="smooth",qt=c(0,1))

cat("Why it works\n")
cat("Calculating principal components...\n")
flush.console()
im.mat <- matrix(images,nr*nc,nf)
im.mat <- im.mat-rowMeans(im.mat)
im.svd <- svd(im.mat)
pcs <- array(im.svd$u, dim=c(nr,nc,nf))
cat("And here they are...\n")
cat("The quadrature signal is assumed to be in the first two slots\n")
flush.console()

if (tolower(.Platform$OS.type) == "windows") {
  windows(width=10, height=2.5*nf)} else {
  x11(width=10, height=2.5*nf)
}
par(mar=rep(0,4))
split.screen(figs=c(nf/2,2))
for (i in 1:nf) {
  screen(i)
  image(pcs[,,i], col=grey256, axes=FALSE, asp=nc/nr, useRaster=TRUE)
}
close.screen(all.screens=TRUE)

#clean up
rm(im.mat, im.svd, pcs)


