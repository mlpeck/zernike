/**********************

Copyright © 2022-2026 Michael Peck <mlpeck54 -at- gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

**************************/


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

# include <cmath>
# include <vector>
# include <RcppArmadillo.h>
# include <boost/multiprecision/float128.hpp>
# include "fastgl_128.h"

using namespace Rcpp;
using namespace arma;
namespace mp = boost::multiprecision;

/*****************
 *
 * elementwise product of two float128 vectors
 *
******************/

std::vector<mp::float128> operator %(const std::vector<mp::float128> x,
                                     const std::vector<mp::float128> y) {
  std::size_t nx = x.size();
  std::vector<mp::float128> prod (nx);
  for (std::size_t i=0; i<nx; ++i) prod[i] = x[i] * y[i];
  return prod;
}

/********************
 *
 * sum the elements of an mp::float128 vector
 *
*********************/

mp::float128 sum_el(std::vector<mp::float128> x) {
  std::size_t nx = x.size();
  mp::float128 sum {0};
  for (std::size_t i=0; i<nx; ++i) sum += x[i];
  return sum;
}



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

void fastgl_xw(const mp::float128& eps,
               std::vector<mp::float128>& xq,
               std::vector<mp::float128>& qwts) {
  
  mp::float128 e2 = eps*eps;
  mp::float128 c1 = (1-e2)/2.;

  std::size_t nq = xq.size();
  
  for (int k=0; k<nq; k++) {
    fastgl::QuadPair p = fastgl::GLPair(nq, k+1);
    xq[k] = c1*(p.x() + 1.0) + e2;
    qwts[k] = c1 * p.weight;
  }
  return;
}
  


mat rzernike_ann_128(const vec& rho, const double& eps, const int& n, const int& m,
                     const std::vector<mp::float128>& xq,
                     const std::vector<mp::float128>& qwts) {

  if (n < m) {
    stop("n < m");
  }
  
  if (2*((n-m)/2) != (n-m)) {
    stop("n,m must be relatively even");
  }
  
  uword nr = rho.n_elem;
  int nz = (n-m)/2 + 1;
  
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
  std::vector<mp::float128> alphal (nz), betal (nz), c (nz);
    
  if (m == 0) {    // know the recursion for this case
    alphal[0] = ak;
    c[0] = e1;
    for (int j = 1; j < nz; j++) {
      alphal[j] = ak;
      betal[j] = j * j * e1l * e1l/4./(2.*j + 1.0)/(2.*j - 1.0);
      c[j] = betal[j] * c[j-1];
    }
  } else { // m > 0
    std::size_t nq = xq.size();
    std::vector<mp::float128> p_k (nq), p_k_plus1 (nq), p_k_minus1 (nq), w (nq);
    mp::float128  inner_product = 0., inner_product_minus1;

    alphal[0] = 0.0;
    for (int i=0; i<nq; i++) {
      w[i] = mp::pow(xq[i], m);
    }
    inner_product = sum_el(w % qwts);
    alphal[0] = sum_el(xq % w % qwts);
    alphal[0] /= inner_product;
    c[0] = inner_product;
    inner_product_minus1 = inner_product;

    for (int i=0; i<nq; i++) {
     p_k[i] = 1.;
     p_k_plus1[i] = xq[i] - alphal[0];
    }

    int k;
    for (k=1; k < (nz-1) ; k++) {
      for (int i=0; i<nq; i++) {
        p_k_minus1[i] = p_k[i];
        p_k[i] = p_k_plus1[i];
      }
      inner_product = sum_el(p_k % p_k % w % qwts);
      alphal[k] = sum_el(xq % p_k % p_k % w % qwts);
      alphal[k] /= inner_product;
      betal[k] = inner_product / inner_product_minus1;
      c[k] = inner_product;
      inner_product_minus1 = inner_product;
      for (int i=0; i<nq; i++) {
        p_k_plus1[i] = (xq[i] - alphal[k]) * p_k[i] - betal[k] * p_k_minus1[i];
      }
    }
    for (int i=0; i<nq; i++) {
      p_k[i] = p_k_plus1[i];
    }
    c[k] = sum_el(p_k % p_k % w % qwts);
  }

  for (int i=0; i<nz; i++) {
    alpha(i) = static_cast<double>(alphal[i]);
    beta(i) = static_cast<double>(betal[i]);
  }

  RZ.col(0).fill(1.0);
  RZ.col(1) = (u - alpha(0));
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
//' @param maxorder the maximum radial polynomial order (defaults to 14).
//' @param nqplus the number of *extra* quadrature nodes (defaults to 6).
//'
//' @return a matrix of Zernike Annular polynomial values evaluated at the input
//'  polar coordinates and all radial orders from
//'  0 through `maxorder` in Fringe sequence, with orthonormal scaling.
//'
//' @details The *radial* polynomials are calculated using recurrence relations
//'  generated numerically using Stieltje's procedure with nominally exact numerical quadrature.
//'  See the documentation for [rzernike_ann()]. A formal presentation is
//'  included in the package documentation.
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
//'     zm <- zpm(rho0, theta0, maxorder=maxorder)
//'     zam <- zapm_128(rhoa, thetaa, eps=eps, maxorder=maxorder)
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
mat zapm_128(const vec& rho, const vec& theta, const double& eps, const int& maxorder=14, const int& nqplus=6) {
  
  int j, k, nmax, mmax = maxorder/2;
  int nz = maxorder/2 + 1;
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
  
  std::size_t nq = nz + nqplus;
  std::vector<mp::float128> xq (nq), qwts (nq);
  mp::float128 epsl(eps);
  fastgl_xw(epsl, xq, qwts);

  //cache values of cos and sin
  
  cosmtheta.col(0) = cos(theta);
  sinmtheta.col(0) = sin(theta);
  for (int m=1; m<mmax; m++) {
    cosmtheta.col(m) = cosmtheta.col(m-1) % cosmtheta.col(0) - sinmtheta.col(m-1) % sinmtheta.col(0);
    sinmtheta.col(m) = sinmtheta.col(m-1) % cosmtheta.col(0) + cosmtheta.col(m-1) % sinmtheta.col(0);
  }
  
  //n=0 zernikes are just the scaled radial zernikes
  
  mat RZ(nr, nz);
  
  RZ = rzernike_ann_128(rho, eps, maxorder, 0, xq, qwts);
  for (int n=0; n<=maxorder; n += 2) {
    k = (n*n)/4 + n;
    zm.col(k) = std::sqrt(n+1.)*RZ.col(n/2);
  }
  
  for (int m=1; m<=mmax; m++) {
    nmax = maxorder - m;
    nz = (nmax - m)/2 + 1;
    mat RZ(nr, nz);
    RZ = rzernike_ann_128(rho, eps, nmax, m, xq, qwts);
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
//' @param maxorder the maximum radial and azimuthal polynomial order (defaults to 14).
//' @param nqplus the number of *extra* quadrature nodes (defaults to 6).
//'
//' @return a matrix of Zernike Annular polynomial values evaluated at the input
//'  polar coordinates and all radial orders from
//'  0 through `maxorder` in ISO/ANSI sequence, with orthonormal scaling.
//'
//' @details The *radial* polynomials are calculated using recurrence relations
//'  generated numerically using Stieltje's procedure with nominally exact numerical quadrature.
//'  See the documentation for [rzernike_ann()]. A formal presentation is
//'  included in the package documentation.
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
//'     zam <- zapm_iso_128(rhoa, thetaa, eps=eps, maxorder=maxorder)
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
mat zapm_iso_128(const vec& rho, const vec& theta, const double& eps, const int& maxorder=14, const int& nqplus=6) {
  
  int j, k, nmax, nz;
  uword nr = rho.size();
  nz = maxorder/2 + 1;
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
  
  int nq = nz + nqplus;
  std::vector<mp::float128> xq (nq), qwts (nq);
  mp::float128 epsl(eps);
  fastgl_xw(epsl, xq, qwts);

  //cache values of cos and sin
  
  cosmtheta.col(0) = cos(theta);
  sinmtheta.col(0) = sin(theta);
  for (int m=1; m<maxorder; m++) {
    cosmtheta.col(m) = cosmtheta.col(m-1) % cosmtheta.col(0) - sinmtheta.col(m-1) % sinmtheta.col(0);
    sinmtheta.col(m) = sinmtheta.col(m-1) % cosmtheta.col(0) + cosmtheta.col(m-1) % sinmtheta.col(0);
  }
  
  //n=0 zernikes are just the scaled radial zernikes
  
  mat RZ(nr, nz);
  
  RZ = rzernike_ann_128(rho, eps, maxorder, 0, xq, qwts);
  for (int n=0; n<=maxorder; n += 2) {
    k = (n*n+2*n)/2;
    zm.col(k) = std::sqrt(n+1.)*RZ.col(n/2);
  }
  
  for (int m=1; m<=maxorder; m++) {
    nmax = maxorder - m % 2;
    nz = (nmax - m)/2 + 1;
    mat RZ(nr, nz);
    RZ = rzernike_ann_128(rho, eps, nmax, m, xq, qwts);
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
