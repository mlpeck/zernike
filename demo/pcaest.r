# demo of PCA fitting
# phase shifted data courtesy of Vladimir Galogaza.

# list of image files
# these are stored in zernike/psidata in the R library, so must expand path name

require(zernike)
	  
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

sa.t <- sconic(diam,roc,lambda=wavelength)
zopt <- get_zoptions()
zopt$plots <- FALSE
zopt$satarget <- sa.t

cat("Performing PSI analysis. Please wait...\n")
flush.console()
system.time(psfit <- psifit(images, phases, psialg="ls", options=zopt))

cat("Zernike fit to wavefront from PSI analysis\n")
flush.console()

plot(psfit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(psfit$wf.smooth),digits=3)))

cat("Performing PCA analysis. Please wait...\n")
system.time(pcfit <- psifit(images, phases, cp=psfit$cp, psialg="pc1", options=zopt))

cat("Wavefront from PCA routine\n")
flush.console()

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(pcfit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(pcfit$wf.smooth),digits=3)))

cat("Performing 'Advanced iterative algorithm'. Please wait...\n")
system.time(aifit <- psifit(images, phases=pcfit$phases, cp=psfit$cp, psialg="aia", options=zopt))

cat("Wavefront from 'AIA'\n")
flush.console()

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(aifit$wf.smooth)
mtext(paste("RMS =",format(pupilrms(aifit$wf.smooth),digits=3)))


# full intercomparison

wfnets <- array(0, dim=c(nr, nc, 3))
wfsmooth <- array(0, dim=c(nr,nc,3))

legends=c("baseline PSI", "PCA", "AIA")

plotn(psfit, pcfit, aifit, labels=legends)
plotn(psfit, pcfit, aifit, labels=legends, wftype="smooth", qt=c(0,1))

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

