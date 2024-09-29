\name{load.images}
\alias{load.images}
\alias{load.pgm}
\title{Read images}
\description{
  Loads image files in jpeg, tiff, png, or raw format.
  \code{load.pgm} provides legacy support for reading
  files in pgm format.
}
\usage{
load.images(files, channels=c(1,0,0), scale=1, FLIP=FALSE)
load.pgm(files, imdiff=NULL)
}
\arguments{
  \item{files}{A vector of character strings with file names}
  \item{channels}{channel weights}
  \item{scale}{scale factor for image resize}
  \item{FLIP}{flip image left for right?}
}
\details{
  set FLIP=TRUE to reverse mirror imaged interferograms.
  
  Any file extension other than jpg, jpeg, tif, tiff, or png is assumed to be in RAW format. Supported raw formats
  are determined by libraw and may not be up to date
}
\value{
  An array containing the contents of the image files.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  \code{load.pgm} is the original \code{load.images}
  included for legacy support of greyscale portable
  anymap files.
}
\keyword{IO}
\keyword{file}
