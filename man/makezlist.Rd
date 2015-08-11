\name{zlist}
\alias{makezlist}
\alias{zlist.fr}
\alias{zmult}
\title{Lists of Zernike polynomial indexes}
\description{
  Ordered lists of Zernike polynomial indexes.
}
\usage{
makezlist(minorder = 2, maxorder = 14)
zlist.fr
zmult(zlist = makezlist())
}
\arguments{
  \item{minorder}{minimum value of \code{n+m}}
  \item{maxorder}{maximum value of \code{n+m}}
  \item{zlist}{a list of the form returned by \code{makezlist}}
}
\details{
  Zernike polynomials are indexed by a radial index \code{n},
  an azimuthal index \code{m}, and include cosine, sine,
  and radial terms. These routines return lists of indexes using
  a popular ordering scheme for Zernike polynomials.
}
\value{
  \code{makezlist} and \code{zlist.fr} return lists with the following components:
  
  \item{n}{radial order}
  \item{m}{azimuthal order}
  \item{t}{one of \code{c("c", "s", "n")}}
  
  \code{zmult} returns a vector the same length as the components
  of \code{zlist}.
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  \code{zlist.fr} is an augmented ``Fringe'' set of Zernike polynomials
     equivalent to makezlist(2,12).
  
  \code{makezlist} returns a complete list of indexes for all orders
  from \code{minorder} through \code{maxorder}, where ``order'' is
  the value of \code{n+m}.
}
\seealso{
  Virtually all high level functions that work with Zernike polynomials
  use these lists. See for example \code{\link{pupil}}, \code{\link{psifit}},
  \code{\link{fftfit}}.
}
\examples{
zlist <- makezlist(2,12)
zcoef <- rnorm(length(zlist))
zcoef		# a vector of normalized Zernike coefficients
zcoef*zmult(zlist)   # Coefficients in conventional representation
sqrt(crossprod(zcoef)) # This is the RMS error of the wavefront 
			# constructed from these Zernikes
}
\keyword{mathematics}
