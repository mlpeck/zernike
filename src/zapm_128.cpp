// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

# include <cmath>
# include <memory>
# include <RcppArmadillo.h>
# include <boost/multiprecision/float128.hpp>
# include "fastgl_128.h"

using namespace Rcpp;
using namespace arma;
namespace mp = boost::multiprecision;

/****************************
 * 
 * nodes and weights from float128 version of fastgl
 * 
 * arguments:
 *   eps obstruction fraction
 *   nq  number of nodes
 *   xq, qwts  nodes and weights for Gauss-Legendre quadrature
 * 
 *   x values are shifted to the interval (eps^2, 1)
 * 
***************************/

int fastgl_xw(const mp::float128& eps, const int& nq, mp::float128 xq[], mp::float128 qwts[]) {
  
  mp::float128 e2 = eps*eps;
  mp::float128 c1 = (1-e2)/2.;
  
  for (int k=0; k<nq; k++) {
    fastgl::QuadPair p = fastgl::GLPair(nq, k+1);
    xq[k] = c1*(p.x() + 1.0) + e2;
    qwts[k] = c1 * p.weight;
  }
  return 0;
}
  


mat rzernike_ann_128(const vec& rho, const double& eps, const int& n, const int& m, const int& nq, const mp::float128 xq[], 
                     const mp::float128 qwts[]) {

  if (n < m) {
    stop("n < m");
  }
  
  if (2*((n-m)/2) != (n-m)) {
    stop("n,m must be relatively even");
  }
  
  uword nr = rho.n_elem;
  int nz = (n-m)/2 + 1;
  int nmax = std::min(2*nz, m+1);
  
  vec u(nr);
  u = rho % rho;

  // return radial Zernikes

  mat RZ(nr, nz);
  vec rm(nr);
  double eps2 = eps*eps;
  double e1 = (1. - eps2);
  mp::float128 epsl(eps);
  mp::float128 e1l(1. - epsl*epsl);
  mp::float128 ak((1. + epsl*epsl)/2.);
  
  rm = pow(rho, m);

  if (nz == 1) {
    RZ.col(0) = std::sqrt(e1/(1-std::pow(eps2, m+1))) * rm;
    return RZ;
  }

  // things we need to calculate for the recurrences
  
  vec alpha(nz), beta(nz);
  mp::float128 nu[2*nz], a[2*nz], b[2*nz];
  mp::float128 alphal[nz], betal[nz], c[nz];
  mp::float128 sigma[3] [2*nz];
  mp::float128 pn[nq], pnp1[nq], pnm1[nq], w[nq];
    
  if (m == 0) {    // know the recursion for this case
    alphal[0] = ak;
    c[0] = e1;
    for (int j = 1; j < nz; j++) {
      alphal[j] = ak;
      betal[j] = j * j * e1l * e1l/4./(2.*j + 1.0)/(2.*j - 1.0);
      c[j] = betal[j] * c[j-1];
    }
  } else { // m > 0
    for (int j = 0; j < 2*nz; j++) {
      a[j] = ak;
      nu[j] = 0.0;
    }
    for (int j = 1; j < 2*nz; j++) {
      b[j] = j * j * e1l * e1l/4./(2.*j + 1.0)/(2.*j - 1.0);
    }
    
    
    // get modified moments
    
    for (int i=0; i<nq; i++) {
      w[i] = pow(xq[i], m);
      nu[0] += w[i] * qwts[i];
      pnm1[i] = 1.0;
      pn[i] = xq[i] - ak;
    }
    
    for (int i=1; i<nmax; i++) {
      for (int j=0; j<nq; j++) {
        nu[i] += pn[j]*w[j]*qwts[j];
      }
      for (int j=0; j<nq; j++) {
        pnp1[j] = (xq[j]-ak) * pn[j] - b[i] * pnm1[j];
        pnm1[j] = pn[j];
        pn[j] = pnp1[j];
      }
    }
    
    //chebyshev algorithm with modified moments
    
    alphal[0] = ak + nu[1]/nu[0];
    betal[0] = nu[0];
    c[0] = nu[0];
    
    for (int l=0; l<(2*nz); l++) {
      sigma[0] [l] = 0.0;
      sigma[1] [l] = nu[l];
    }
  
    for (int k=1; k<nz; k++) {
      for (int l=k; l<(2*nz-k); l++) {
        sigma[2][l] = sigma[1][l+1] + (ak - alphal[k-1]) * sigma[1][l] - 
                        betal[k-1] * sigma[0][l] + b[l] * sigma[1][l-1];
      }
      alphal[k] = ak + sigma[2][k+1]/sigma[2][k] - sigma[1][k]/sigma[1][k-1];
      betal[k] = sigma[2][k]/sigma[1][k-1];
      c[k] = betal[k] * c[k-1];
      for (int l=k; l<(2*nz-k); l++) {
        sigma[0][l] = sigma[1][l];
        sigma[1][l] = sigma[2][l];
      }
    }
  }
  for (int i=0; i<nz; i++) {
    alpha(i) = static_cast<double>(alphal[i]);
    beta(i) = static_cast<double>(betal[i]);
  }

  RZ.col(0).fill(1.0);
  RZ.col(1) = (u - alpha(0)) % RZ.col(0);
  for (int i = 1; i < (nz-1); i++) {
    RZ.col(i+1) = (u - alpha(i)) % RZ.col(i) - beta(i) * RZ.col(i-1);
  }
  for (int i=0; i<nz; i++) {
    RZ.col(i) = RZ.col(i) % rm * std::sqrt(e1/(2. * i + m + 1.)/static_cast<double>(c[i]));
  }
  
  return RZ;
}


/*****************
 * 
 * Create matrix of Zernike Annular polynomials
 * in extended Fringe index scheme.
 * 
 * 
******************/

//' Zernike Annular polynomials, extended precision version
//'
//' Create a matrix of Zernike Annular polynomial values in
//' extended Fringe sequence for a set of polar coordinates.
//'
//' @param rho a vector of radial coordinates with eps <= rho <= 1.
//' @param theta a vector of angular coordinates, in radians.
//' @param eps the obstruction fraction 0 <= eps < 1.
//' @param maxorder the maximum radial polynomial order (defaults to 12).
//' @param nq the number of quadrature points for numerical integration
//'
//' @return a matrix of Zernike Annular polynomial values evaluated at the input
//'  polar coordinates and all radial orders from
//'  0 through `maxorder` in Fringe sequence, with orthonormal scaling.
//'
//' @details The *radial* polynomials are calculated using recurrence relations
//'  generated numerically using chebyshev's algorithm with modified moments.
//'  See the documentation for [rzernike_ann()]. A formal presentation will be
//'  included in a future release.
//' @examples
//'   sample_az_128 <- function(maxorder=12, eps=0.33, col=rev(zernike::rygcb(400)), addContours=TRUE, cscale=TRUE) {
//'   
//'     ## get coordinates for unobstructed and obstructed apertures
//'     cpa <- cp.default
//'     cpa$obstruct <- eps
//'     prt <- pupil.rhotheta(nrow.default,ncol.default,cp.default)
//'     prta <- pupil.rhotheta(nrow.default,ncol.default,cp=cpa)
//'     rho0 <- prt$rho[!is.na(prt$rho)]
//'     theta0 <- prt$theta[!is.na(prt$theta)]
//'     rhoa <- prta$rho[!is.na(prta$rho)]
//'     thetaa <- prta$theta[!is.na(prta$theta)]
//'     
//'     ## fill up matrixes of Zernikes and Annular Zernikes
//'     
//'     zm <- zpmC(rho0, theta0, maxorder=maxorder)
//'     zam <- zapm_128(rhoa, thetaa, eps=eps, maxorder=maxorder, nq=maxorder/2+5)
//'     
//'     ## pick a column at random and look up its index pair
//'     
//'     zlist <- makezlist(0, maxorder)
//'     i <- sample(2:ncol(zm), 1)
//'     n <- zlist$n[i]
//'     m <- zlist$m[i]
//'     
//'     ## fill up the wavefront representations and plot them
//'     
//'     wf0 <- prt$rho
//'     wf0[!is.na(wf0)] <- zm[,i]
//'     class(wf0) <- "pupil"
//'     
//'     wfa <- prta$rho
//'     wfa[!is.na(wfa)] <- zam[,i]
//'     class(wfa) <- "pupil"
//'     
//'     plot(wf0, cp=cp.default, col=col, addContours=addContours, cscale=cscale)
//'     mtext(paste("Zernike, n =", n, " m =", m))
//'     
//'     x11()
//'     plot(wfa, cp=cpa, col=col, addContours=addContours, cscale=cscale)
//'     mtext(paste("Annular Zernike, n =", n, " m =", m))
//'     
//'     ## return Zernike matrices and wavefronts invisibly
//'     ## just in case user wants to do something with them
//'     
//'     invisible(list(zm=zm, wf0=wf0, zam=zam, wfa=wfa))
//'   }
//'
//'   sample_az_128()
//'
//' @md
// [[Rcpp::export]]
mat zapm_128(const vec& rho, const vec& theta, const double& eps, const int& maxorder=12, const int& nq=21) {
  
  int j, k, n0, nmax, nz, mmax = maxorder/2;
  uword nr = rho.size();
  int ncol = (mmax+1)*(mmax+1);
  mat cosmtheta(nr, mmax), sinmtheta(nr, mmax);
  mat zm(nr, ncol);
  
    //do some rudimentary error checking
  
  if (rho.size() != theta.size()) {
    stop("Numeric vectors must be same length");
  }
  if ((maxorder % 2) != 0) {
    stop("maxorder must be even");
  }
  
  //good enough
  
  // get points and weights for quadrature
  
  mp::float128 xq[nq], qwts[nq];
  mp::float128 epsl(eps);
  fastgl_xw(epsl, nq, xq, qwts);

  //cache values of cos and sin
  
  cosmtheta.col(0) = cos(theta);
  sinmtheta.col(0) = sin(theta);
  for (int m=1; m<mmax; m++) {
    cosmtheta.col(m) = cosmtheta.col(m-1) % cosmtheta.col(0) - sinmtheta.col(m-1) % sinmtheta.col(0);
    sinmtheta.col(m) = sinmtheta.col(m-1) % cosmtheta.col(0) + cosmtheta.col(m-1) % sinmtheta.col(0);
  }
  
  //n=0 zernikes are just the scaled radial zernikes
  
  nz = maxorder/2 + 1;
  mat RZ(nr, nz);
  
  RZ = rzernike_ann_128(rho, eps, maxorder, 0, nq, xq, qwts);
  for (int n=0; n<=maxorder; n += 2) {
    k = (n*n)/4 + n;
    zm.col(k) = std::sqrt(n+1.)*RZ.col(n/2);
  }
  
  for (int m=1; m<=mmax; m++) {
    nmax = maxorder - m;
    nz = (nmax - m)/2 + 1;
    mat RZ(nr, nz);
    RZ = rzernike_ann_128(rho, eps, nmax, m, nq, xq, qwts);
    j = 0;
    for (int n=m; n<= nmax; n += 2) {
      k = ((n+m)*(n+m))/4 + n - m;
      zm.col(k) = std::sqrt(2.*(n+1.)) * RZ.col(j) % cosmtheta.col(m-1);
      k++;
      zm.col(k) = std::sqrt(2.*(n+1.)) * RZ.col(j) % sinmtheta.col(m-1);
      j++;
    }
  }
    
  return zm;
}

/*****************
 * 
 * Create matrix of Zernike Annular polynomials
 * in ISO/ANSI index scheme.
 * 
 * 
******************/

//' Zernike Annular polynomials, ISO ordering - extended precision version
//'
//' Create a matrix of Zernike Annular polynomial values in
//' ISO/ANSI sequence for a set of polar coordinates.
//'
//' @param rho a vector of radial coordinates with eps <= rho <= 1.
//' @param theta a vector of angular coordinates, in radians.
//' @param eps the obstruction fraction 0 <= eps < 1.
//' @param maxorder the maximum radial and azimuthal polynomial order (defaults to 12).
//' @param nq the number of quadrature points for numerical integration
//'
//' @return a matrix of Zernike Annular polynomial values evaluated at the input
//'  polar coordinates and all radial orders from
//'  0 through `maxorder` in ISO/ANSI sequence, with orthonormal scaling.
//'
//' @details The *radial* polynomials are calculated using recurrence relations
//'  generated numerically using chebyshev's algorithm with modified moments.
//'  See the documentation for [rzernike_ann()]. A formal presentation will be
//'  included in a future release.
//'
//' @examples
//'   sample_az_iso_128 <- function(maxorder=12, eps=0.33, col=rev(zernike::rygcb(400)), addContours=TRUE, cscale=TRUE) {
//'   
//'     ## get coordinates for unobstructed and obstructed apertures
//'     cpa <- cp.default
//'     cpa$obstruct <- eps
//'     prt <- pupil.rhotheta(nrow.default,ncol.default,cp.default)
//'     prta <- pupil.rhotheta(nrow.default,ncol.default,cp=cpa)
//'     rho0 <- prt$rho[!is.na(prt$rho)]
//'     theta0 <- prt$theta[!is.na(prt$theta)]
//'     rhoa <- prta$rho[!is.na(prta$rho)]
//'     thetaa <- prta$theta[!is.na(prta$theta)]
//'     
//'     ## fill up matrixes of Zernikes and Annular Zernikes
//'     
//'     zm <- zpm_cart(x=rho0*cos(theta0), y=rho0*sin(theta0), maxorder=maxorder)
//'     zam <- zapm_iso_128(rhoa, thetaa, eps=eps, maxorder=maxorder, nq=maxorder/2+5)
//'     
//'     ## pick a column at random and look up its index pair
//'     
//'     zlist <- makezlist.iso(maxorder)
//'     i <- sample(2:ncol(zm), 1)
//'     n <- zlist$n[i]
//'     m <- zlist$m[i]
//'     
//'     ## fill up the wavefront representations and plot them
//'     
//'     wf0 <- prt$rho
//'     wf0[!is.na(wf0)] <- zm[,i]
//'     class(wf0) <- "pupil"
//'     
//'     wfa <- prta$rho
//'     wfa[!is.na(wfa)] <- zam[,i]
//'     class(wfa) <- "pupil"
//'     
//'     plot(wf0, cp=cp.default, col=col, addContours=addContours, cscale=cscale)
//'     mtext(paste("Zernike, n =", n, " m =", m))
//'     
//'     x11()
//'     plot(wfa, cp=cpa, col=col, addContours=addContours, cscale=cscale)
//'     mtext(paste("Annular Zernike, n =", n, " m =", m))
//'     
//'     ## return Zernike matrices and wavefronts invisibly
//'     ## just in case user wants to do something with them
//'     
//'     invisible(list(zm=zm, wf0=wf0, zam=zam, wfa=wfa))
//'   }
//'
//'   sample_az_iso_128()
//'
//' @md
// [[Rcpp::export]]
mat zapm_iso_128(const vec& rho, const vec& theta, const double eps, const int& maxorder=12, const int& nq=21) {
  
  int j, k, nmax, nz;
  uword nr = rho.size();
  int ncol = (maxorder+1)*(maxorder+2)/2;
  mat cosmtheta(nr, maxorder), sinmtheta(nr, maxorder);
  mat zm(nr, ncol);
  
    //do some rudimentary error checking
  
  if (rho.size() != theta.size()) {
    stop("Numeric vectors must be same length");
  }
  if (maxorder < 1) {
    stop("maxorder must be >= 1");
  }
  //good enough
  
  // get points and weights for quadrature
  
  mp::float128 xq[nq], qwts[nq];
  mp::float128 epsl(eps);
  fastgl_xw(epsl, nq, xq, qwts);

  //cache values of cos and sin
  
  cosmtheta.col(0) = cos(theta);
  sinmtheta.col(0) = sin(theta);
  for (int m=1; m<maxorder; m++) {
    cosmtheta.col(m) = cosmtheta.col(m-1) % cosmtheta.col(0) - sinmtheta.col(m-1) % sinmtheta.col(0);
    sinmtheta.col(m) = sinmtheta.col(m-1) % cosmtheta.col(0) + cosmtheta.col(m-1) % sinmtheta.col(0);
  }
  
  //n=0 zernikes are just the scaled radial zernikes
  
  nz = maxorder/2 + 1;
  mat RZ(nr, nz);
  
  RZ = rzernike_ann_128(rho, eps, maxorder, 0, nq, xq, qwts);
  for (int n=0; n<=maxorder; n += 2) {
    k = (n*n+2*n+2)/2 - 1;
    zm.col(k) = std::sqrt(n+1.)*RZ.col(n/2);
  }
  
  for (int m=1; m<=maxorder; m++) {
    nmax = maxorder - m % 2;
    nz = (nmax - m)/2 + 1;
    mat RZ(nr, nz);
    RZ = rzernike_ann_128(rho, eps, nmax, m, nq, xq, qwts);
    j = 0;
    for (int n=m; n<= nmax; n += 2) {
      k = (n*n + 2*n - m)/2;
      zm.col(k) = std::sqrt(2.*(n+1.)) * RZ.col(j) % sinmtheta.col(m-1);
      k += m;
      zm.col(k) = std::sqrt(2.*(n+1.)) * RZ.col(j) % cosmtheta.col(m-1);
      j++;
    }
  }
    
  return zm;
}
