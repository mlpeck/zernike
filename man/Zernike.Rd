\name{Zernike}
\alias{Zernike}
\alias{rzernike}
\alias{drzernike}
\title{Zernike Polynomials}
\description{
  Routines for creating and manipulating Zernike polynomials.
}
\usage{
Zernike(rho, theta, n, m, t)
rzernike(rho, n, m)
drzernike(rho, n, m)
}
\arguments{
  \item{rho}{normalized radius, \eqn{0 <= rho <= 1}}
  \item{theta}{angular coordinate}
  \item{n}{radial polynomial order}
  \item{m}{azimuthal order}
  \item{t}{character for trig function: one of c("n", "c", "s")}
}
\references{
  Born, M. and Wolf, E. 1999, \emph{Principles of Optics, 7th Edition},
  Cambridge University Press, chapter 9 and appendix VII.

  Noll, R.J. 1976, \bold{Zernike polynomials and atmospheric turbulence}, \emph{J. Opt. Soc. Am.},
  Vol. 66, No. 3, p. 207.

  \url{http://wyant.opt-sci.arizona.edu/zernikes/zernikes.htm}

  \url{http://mathworld.wolfram.com/ZernikePolynomial.html}
}
\author{M.L. Peck \email{mpeck1@ix.netcom.com}}
\note{
  These functions return Zernikes scaled such that they form an orthonormal basis set
  for the space of functions defined on the unit circle. Note that this is not the most
  commonly used definition (as given e.g. in \emph{Born and Wolf}). The definition I use is
  often associated with \emph{Noll (1976)}.

  The function \code{\link{zmult}} can be used to convert between normalized
  and conventionally defined vectors of Zernike coefficients.
  
  The basic low level functions \code{rzernike} and \code{drzernike}
  use numerically stable recurrence relationships for the radial Zernikes.
}
\seealso{
  \code{\link{makezlist}},
  \code{\link{zlist.fr}},
  \code{\link{zmult}},
  \code{\link{zpm}},
  \code{\link{pupil}},
  \code{\link{pupilrms}},
  \code{\link{pupilpv}},
  \code{\link{strehlratio}}.
}
\examples{
Zernike(1, 0, 4, 0, "n")	# == sqrt(5)

# A slightly more complex example

rho <- seq(0, 1, length = 101)
theta <- rep(0, 101)

plot(rho, Zernike(rho, theta, 6, 0, "n"), type="l", 
  ylim=c(-3.5,3.5), main="Some 6th order Zernike Polynomials")
lines(rho, Zernike(rho, theta, 5, 1, "c"), lty=2)
lines(rho, Zernike(rho, theta, 4, 2, "c"), lty=3)
lines(rho, Zernike(rho, theta, 3, 3, "c"), lty=4)
}
\keyword{mathematics}
