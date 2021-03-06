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
  nrow.default <- 640
  ncol.default <- nrow.default
  cp.default <- list(xc=320.5, yc=320.5, rx=319.5, ry=319.5, obstruct=0)

  pupil(zcoef=NULL, maxorder=14L, isoseq=FALSE,
      phi=0, piston=NULL,
      nrow=nrow.default, ncol=ncol.default, 
      cp=cp.default)
  pupil.arb(zcoef=NULL, zlist=makezlist(), 
          phi=0, piston=0,
          nrow=nrow.default, ncol=ncol.default, 
          cp=cp.default)
}
\arguments{
  \item{zcoef}{A vector of Zernike coefficients}
  \item{maxorder}{Maximum Zernike polynomial order}
  \item{zlist}{List of indexes the same length as zcoef}
  \item{isoseq}{ZPs in ISO/ANSI sequence}
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

  \code{pupil.arb} will accept an arbitrary list of Zernike indexes.

  \code{pupil} requires a complete set of Zernikes as returned
  by \code{\link{zpm}} or \code{\link{zpm_cart}}.
  
  Default values \code{nrow.default}, \code{ncol.default} and \code{cp.default}
  are now (as of version 3.7.0) stored in the package environment. These
  can be used or overridden by setting values in the user`s global environment.
  
  Version 3.7.0 of this package makes some fairly significant changes in \code{pupil} from previous releases.
  First, zlist is no longer used to indicate the set of Zernikes used. Instead the 
  maximum radial order is passed as \code{maxorder}, defaulting to 14.
  The complete set of Zernike polynomial values including piston is
  calculated with either \code{\link{zpm}} or \code{\link{zpm_cart}} if
  \code{isoseq} is \code{TRUE}. The vector \code{zcoef}, if non-null,
  must be either the same length as the number of columns in the zernike matrix with
  a null value for \code{piston}, or one shorter with a non-null value
  for \code{piston}. The vector \code{zcoef.net} returned by
  \code{\link{wf_net}} has the first, piston, element removed so a
  non-null piston value must be passed as an argument.
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
  \code{\link{zpm}},
  \code{\link{zpm_cart}},
  \code{\link{wf_net}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}},
  \code{\link{strehlratio}},
  \code{\link{pupil.pars}},
  \code{\link{circle.pars}}.
}
\examples{
wf <- pupil(zcoef=rnorm(length(makezlist(0, 14)$n), 0, 0.01))
plot(wf, addContours=FALSE)
summary(wf)
}
\keyword{Graphics}
