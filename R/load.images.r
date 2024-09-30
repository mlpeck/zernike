load.images <- function(files, channels=c(1,0,0), scale=1, FLIP=FALSE) {
  channels <- channels/sum(channels)
  ext <- strsplit(files[1], ".", fixed=TRUE)[[1]]
  ext <- tolower(ext[length(ext)])
  switch(ext, 
    jpg = , 
    jpeg = {
      if (!require(jpeg)) {
        stop("Install package jpeg to read jpeg files")
      } else {
        fn <- readJPEG
        fromraw <- FALSE
      }
    },
    tif = , 
    tiff = {
      if (!require(tiff)) {
        stop("Install package tiff to read tiff files")
      } else {
        fn <- readTIFF
        fromraw <- FALSE
      }
    },
    png = {
      if (!require(png)) {
        stop("Install package png to read png files")
      } else {
        fn <- readPNG
        fromraw <- FALSE
      }
    },
    dng = ,
    nef = {
      fn <- readraw 
      fromraw <- TRUE
    }, ## if no match assume raw
    {
      fn <- readraw 
      fromraw <- TRUE
    }
  )
    
  readf <- function(fn, file, channels, fromraw=FALSE) {
    if (fromraw) {
      im <- fn(file, channels)
    } else {
      im <- fn(file)
      dims <- dim(im)
      if (length(dims) > 2) {
        im0 <- matrix(0, dims[1], dims[2])
        for (i in seq_along(channels)) {
          im0 <- im0 + im[,,i]*channels[i]
        }
        im <- im0
      }
    }
    if (scale != 1) {
      im <- rescale(im, scale)
    }
    if (!fromraw) {
      t(im[nrow(im):1,])
    } else {
      im[, ncol(im):1]
    }
  }
    
  if (length(files)==1) {
    images <- readf(fn, files, channels, fromraw)
    if(FLIP) {
      images <- images[nrow(images):1, ]
    }
  } else {
    im1 <- readf(fn, files[1], channels, fromraw)
    nr <- nrow(im1)
    nc <- ncol(im1)
    images <- array(0, dim=c(nr, nc, length(files)))
    images[,,1] <- im1
    for (i in 2:length(files)) {
      images[,,i] <- readf(fn, files[i], channels, fromraw)
    }
    if(FLIP) {
      images <- images[nr:1,,]
    }
  }
  ri <- range(images)
  images <- (images - ri[1])/(ri[2]-ri[1])
  images
}
