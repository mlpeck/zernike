\name{pupil}
\alias{pupil}
\alias{pupil.arb}
\title{Pupils and wavefronts}
\description{
  Create a pupil object and optionally fill it with a wavefront.
  For our purposes a ``pupil'' is defined to be a matrix representation
  of a circular or annular aperture. Simple plot and summary methods
  are also provided.
}
\usage{
  pupil(zcoef=NULL, maxorder=12L, isoseq=FALSE,
      phi=0, piston=0,
      nrow=640, ncol=nrow, 
      cp=list(xc=320.5,yc=320.5,rx=319.5,ry=319.5,obstruct=0))
  pupil.arb(zcoef=NULL, zlist=makezlist(), 
          phi=0, piston=0,
          nrow=640, ncol=nrow, 
          Cp=list(xc=320.5,yc=320.5,rx=319.5,ry=319.5,obstruct=0))
}
\arguments{
  \item{zcoef}{A vector of Zernike coefficients}
  \item{maxorder}{Maximum Zernike polynomial order}
  \item{isoseq}{ZPs in ISO/ANSI sequence}
  \item{zlist}{List of indexes the same length as \code{zcoef}}
  \item{phi}{Amount to rotate image, in degrees}
  \item{piston}{Constant to add to wavefront values}
  \item{nrow}{Number of rows in output matrix}
  \item{ncol}{Number of columns in output matrix}
  \item{cp}{A list with items
    \code{xc} - x coordinate of central pixel,
    \code{yc} - y coordinate of central pixel,
    \code{rx} - x radius in pixels,
    \code{ry} - y radius in pixels,
    \code{obstruct} - central obstruction fraction.
  }
}
\details{
  \code{\link{plot.pupil}} and \code{\link{summary.pupil}} provide
  simple plot and summary methods for objects of class \code{"pupil"}.

  \code{pupil.arb} will accept an arbitrary list of Zernikes.

  \code{pupil} requires a complete set of Zernikes as returned
  by \code{\link{makezlist}}.
}
\value{
  A matrix of size \code{nrow} x \code{ncol}.
  The matrix is assigned to the class \code{"pupil"}.
  \code{NA}s are used to fill the matrix outside the defined
  area of the pupil.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  The parameter \code{cp} is used to define the dimensions of the pupil.
  See \code{\link{pupil.pars}} for details.
  
}

\seealso{
  \code{\link{Zernike}},
  \code{\link{makezlist}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}},
  \code{\link{strehlratio}},
  \code{\link{pupil.pars}},
  \code{\link{circle.pars}}.
}
\examples{
wf <- pupil(zcoef=rnorm(length(makezlist()$n), 0, 0.01))
plot(wf, addContours=FALSE)
summary(wf)
}
\keyword{Graphics}
