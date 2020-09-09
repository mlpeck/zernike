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


## load image files (jpeg or tiff format) into an array

load.images <- function(files, channels=c(1,0,0), scale=1, FLIP=FALSE) {
    ext <- strsplit(files[1], ".", fixed=TRUE)[[1]]
    ext <- tolower(ext[length(ext)])
    switch(ext, jpg = , jpeg = {fn = readjpeg; fromraw=FALSE}, 
      tif = , tiff = {fn=readtiff; fromraw=FALSE},
      {fn = readraw; fromraw=TRUE})

    readf <- function(fn, file, channels, fromraw=FALSE) {
        im <- fn(file, channels)
        if (scale != 1) im <- rescale(im, scale)
        if (!fromraw) {
          t(im[nrow(im):1,])
        } else {
          im[, ncol(im):1]
        }
    }
    if (length(files)==1) {
        images <- readf(fn, files, channels, fromraw)
        ri <- range(images)
        images <- (images-ri[1])/(ri[2]-ri[1])
        if(FLIP) images <- images[nrow(images):1, ]
        return(images)
    }
    im1 <- readf(fn, files[1], channels, fromraw)
    nr <- nrow(im1)
    nc <- ncol(im1)
    images <- array(0, dim=c(nr, nc, length(files)))
    images[,,1] <- im1
    for (i in 2:length(files)) {
        images[,,i] <- readf(fn, files[i], channels, fromraw)
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
  res <- .C(resize_image, as.double(im), 
			as.integer(ncol(im)), as.integer(nrow(im)),
			im.out = double(nr.new*nc.new),
			as.integer(nc.new), as.integer(nr.new),
			ret = integer(1), NAOK=TRUE)
  im.out <- matrix(res$im.out, nr.new, nc.new)
  im.out
}
  
## read a jpeg file

readjpeg <- function(filename, channels){
    res <- .C(read_jpg_img_info, as.character(filename),
	          width=integer(1), height=integer(1), depth=integer(1),
	          ret=integer(1))
    if (res$ret < 0) stop(if (res$ret == -1) "Cannot open file." else "Internal error")
    width <- res$width
    height <- res$height
    depth <- res$depth
    res <- .C(read_jpg_img, as.character(filename),
                  as.integer(width), as.integer(height), as.integer(depth),
                  as.double(channels),
	          image=double(res$width * res$height),
	          ret=integer(1))
    image <- matrix(res$image, height, width)
    image
}

## read a tiff file
readtiff <- function(filename, channels){
    res <- .C(read_tiff_img_info, as.character(filename),
	          width=integer(1), height=integer(1), depth=integer(1),
	          ret=integer(1))
    if (res$ret < 0) stop(if (res$ret == -1) "Cannot open file." else "Internal error")
    width <- res$width
    height <- res$height
    depth <- res$depth
    res <- .C(read_tiff_img, as.character(filename),
                  as.integer(width), as.integer(height), as.integer(depth),
                  as.double(channels),
	          image=double(res$width * res$height),
	          ret=integer(1))
    image <- matrix(res$image, height, width)
    image
}

