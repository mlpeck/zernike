% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{zapm_iso}
\alias{zapm_iso}
\title{Zernike Annular polynomials, ISO ordering}
\usage{
zapm_iso(rho, theta, eps, maxorder = 12L, nq = maxorder/2L + 5L)
}
\arguments{
\item{rho}{a vector of radial coordinates with eps <= rho <= 1.}

\item{theta}{a vector of angular coordinates, in radians.}

\item{eps}{the obstruction fraction 0 <= eps < 1.}

\item{maxorder}{the maximum radial and azimuthal polynomial order (defaults to 12).}

\item{nq}{the number of quadrature points for numerical integration}
}
\value{
a matrix of Zernike Annular polynomial values evaluated at the input
polar coordinates and all radial orders from
0 through \code{maxorder} in ISO/ANSI sequence, with orthonormal scaling.
}
\description{
Create a matrix of Zernike Annular polynomial values in
ISO/ANSI sequence for a set of polar coordinates.
}
\details{
The \emph{radial} polynomials are calculated using recurrence relations
generated numerically using chebyshev's algorithm with modified moments.
See the documentation for \code{\link[=rzernike_ann]{rzernike_ann()}}. A formal presentation is
included in the package documentation.
}
\examples{
  sample_az_iso <- function(maxorder=12, eps=0.33, col=rev(zernike::rygcb(400)), addContours=TRUE, cscale=TRUE) {
  
    ## get coordinates for unobstructed and obstructed apertures
    cpa <- cp.default
    cpa$obstruct <- eps
    prt <- pupil.rhotheta(nrow.default,ncol.default,cp.default)
    prta <- pupil.rhotheta(nrow.default,ncol.default,cp=cpa)
    rho0 <- prt$rho[!is.na(prt$rho)]
    theta0 <- prt$theta[!is.na(prt$theta)]
    rhoa <- prta$rho[!is.na(prta$rho)]
    thetaa <- prta$theta[!is.na(prta$theta)]
    
    ## fill up matrixes of Zernikes and Annular Zernikes
    
    zm <- zpm_cart(x=rho0*cos(theta0), y=rho0*sin(theta0), maxorder=maxorder)
    zam <- zapm_iso(rhoa, thetaa, eps=eps, maxorder=maxorder, nq=maxorder/2+5)
    
    ## pick a column at random and look up its index pair
    
    zlist <- makezlist.iso(maxorder)
    i <- sample(2:ncol(zm), 1)
    n <- zlist$n[i]
    m <- zlist$m[i]
    
    ## fill up the wavefront representations and plot them
    
    wf0 <- prt$rho
    wf0[!is.na(wf0)] <- zm[,i]
    class(wf0) <- "pupil"
    
    wfa <- prta$rho
    wfa[!is.na(wfa)] <- zam[,i]
    class(wfa) <- "pupil"
    
    plot(wf0, cp=cp.default, col=col, addContours=addContours, cscale=cscale)
    mtext(paste("Zernike, n =", n, " m =", m))
    
    x11()
    plot(wfa, cp=cpa, col=col, addContours=addContours, cscale=cscale)
    mtext(paste("Annular Zernike, n =", n, " m =", m))
    
    ## return Zernike matrices and wavefronts invisibly
    ## just in case user wants to do something with them
    
    invisible(list(zm=zm, wf0=wf0, zam=zam, wfa=wfa))
  }

  sample_az_iso()

}
