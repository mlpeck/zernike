## read an image in 'portable any map format'. Included for legacy support


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
    images
}

