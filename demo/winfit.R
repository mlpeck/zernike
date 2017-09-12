# moving window fits

require(zernike)

# check to see if object "images" is in workspace. If not, read in data.

if (!exists("images"))
{
fpath <- file.path(find.package(package="zernike"), "psidata")
files <- scan(file.path(fpath, "files.txt"), what="character")
for (i in 1:length(files)) files[i] <- file.path(fpath, files[i])

# load the images into an array

images <- load.images(files)

# parameters for this run

source(file.path(fpath, "parameters.txt"))

# phase shifts

phases <- wrap((0:(frames.per.cycle-1))/frames.per.cycle*2*pi)
phases <- switch(ps.dir, ccw = -phases, cw = phases, phases)

# target SA coefficients for numerical null.

sa.t <- sconic(diam,roc,lambda=wavelength)

}


n.w <- dim(images)[3]-frames.per.cycle+1
zlist <- makezlist(2,14)
zcoefs.w <- NULL
cp.w <- NULL
for (i in 1:n.w) {
  # fit each cycle
	if (tolower(.Platform$OS.type) == "windows") windows() else x11()
    temp <- psifit(images[,,i:(i+frames.per.cycle-1)],
      phases[1:frames.per.cycle],satarget=sa.t,zlist=zlist)
  # this assigns the fit to an object named psfit.1, etc.
    assign(paste("psfit", i, sep="."), temp)
  # get some information about the fit and store it in a matrix
    zcoefs.w <- rbind(zcoefs.w, c(sqrt(crossprod(temp$zcoef.net)), temp$fit[-1]))
    cp.w <- rbind(cp.w, c(temp$cp$xc, temp$cp$yc, temp$cp$rx))
}
# make these data frames
zcoefs.w <- data.frame(zcoefs.w)
colnames(zcoefs.w) <- c("rms", "tilt.x", "tilt.y", "defocus", "astig.x", "astig.y",
        "coma.x","coma.y","SA","trefoil.x","trefoil.y","astig5.x","astig5.y",
        "coma5.x","coma5.y","SA5",paste("Z",16:63,sep=""))

cp.w <- data.frame(cp.w)
colnames(cp.w) <- c("x.c", "y.c", "rxy")

# print out mean and standard deviation of cp.w

colMeans(cp.w)
sapply(cp.w, sd)

#subtract asphere from Zernike coefficients

zcoefs.w[,9] <- zcoefs.w[,9]-sa.t[1]
zcoefs.w[,16] <- zcoefs.w[,16]-sa.t[2]

zorder <- (zlist$n+zlist$m)[-(1:3)]

if (tolower(.Platform$OS.type) == "windows") windows() else x11()
plot(zorder, zcoefs.w[1,5:64], xlab="Zernike order", ylab="Coefficient value",
  type="n", ylim=range(zcoefs.w[,5:64]))
for (i in 1:nrow(zcoefs.w))
  points(zorder, zcoefs.w[i,5:64], pch=zlist$m[-(1:3)]+1, col=zlist$m[-(1:3)]+1)
  legend(12, 0.085, legend=0:max(zlist$m), pch=1:(max(zlist$m)+1),
    col=1:(max(zlist$m)+1), title="Azimuthal order")
title(main="Zernike coefficients")

rm(temp)
