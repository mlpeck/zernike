//
//	Copyright (C) 2019 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// Fill a matrix with Zernike polynomial values

// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

NumericMatrix zpmC(const NumericVector& rho, const NumericVector& theta, const int& maxorder = 12) {

  int m, n, n0, mmax = maxorder/2;
  unsigned int i, imax = rho.size();
  int order, nm, nm1mm1, nm1mp1, nm2m;
  int ncol = (mmax+1)*(mmax+1);
  double a0;
  double cosmtheta[mmax], sinmtheta[mmax];
  NumericMatrix zm(imax, ncol);
  
  //do some rudimentary error checking
  
  if (rho.size() != theta.size()) stop("Numeric vectors must be same length");
  if ((maxorder % 2) != 0) stop("maxorder must be even");
  
  //good enough
  
  
  for (i=0; i<imax; i++) {
    
    //cache values of cos and sin
    cosmtheta[0] = std::cos(theta[i]);
    sinmtheta[0] = std::sin(theta[i]);
    for (m=1; m<mmax; m++) {
      cosmtheta[m] = cosmtheta[m-1]*cosmtheta[0] - sinmtheta[m-1]*sinmtheta[0];
      sinmtheta[m] = sinmtheta[m-1]*cosmtheta[0] + cosmtheta[m-1]*sinmtheta[0];
    }
    
    zm(i, 0) = 1.0;                     //piston term
    zm(i, 3) = 2. * rho[i] * rho[i] - 1.; //defocus
    
    // now fill in columns with m=n for n>0
    
    for (m=1; m <= mmax; m++) {
      zm(i, m*m) = rho[i] * zm(i, (m-1)*(m-1));
    }

    // non-symmetric terms
    
    for (order=4; order<=maxorder; order+=2) {
      for (m=order/2-1; m>0; m--) {
        n=order-m;
        nm = order*order/4 + n - m;
        nm1mm1 = (order-2)*(order-2)/4 + n - m;
        nm1mp1 = nm - 2;
        nm2m = nm1mm1 - 2;
        zm(i, nm) = rho[i]*(zm(i, nm1mm1) + zm(i, nm1mp1)) - zm(i, nm2m);
      }
      
      // m=0 (symmetric) term
      nm = order*order/4 + order;
      nm1mp1 = nm-2;
      nm2m = (order-2)*(order-2)/4+order-2;
      zm(i, nm) = 2.*rho[i]*zm(i, nm1mp1) - zm(i, nm2m);
    }

    // now multiply each column by normalizing factor and cos, sin
    
    n0 = 1;
    for (order=2; order <= maxorder; order+=2) {
      for(m=order/2; m>0; m--) {
        n=order-m;
        a0 = sqrt(2.*(n+1));
        zm(i, n0+1) = a0*sinmtheta[m-1]*zm(i, n0);
        zm(i, n0) *= a0*cosmtheta[m-1];
        n0 += 2;
      }
      n = order;
      zm(i, n0) *= sqrt(n+1.);
      n0++;
    }
  }
  return zm;
}

/*****************
 * 
 * Approximate Zernike Annular polynomials by numerically orthogonalizing
 * sub-matrixes for each m separately using thin QR decomposition.
 * This version is for extended Fringe set order and polar coordinates.
 * 
******************/

// [[Rcpp::export]]
mat zapmC(const NumericVector& rho, const NumericVector& theta, const int& maxorder=12) {
  
  uword nrow = rho.size();
  int mmax = maxorder/2;
  int ncol = (mmax+1)*(mmax+1);
  int i, j, m, nj;
  double zpnorm = std::sqrt((double) nrow);
  mat zm(nrow, ncol), annzm(nrow, ncol);
  
  //do some rudimentary error checking
  
  if (rho.size() != theta.size()) stop("Numeric vectors must be same length");
  if ((maxorder % 2) != 0) stop("maxorder must be even");
  
  //good enough
  
  zm = as<arma::mat>(zpmC(rho, theta, maxorder));
  
  // for each azimuthal index m find the column indexes
  // of the zp matrix having that value of m. That
  // subset is what we need to orthogonalize.
  // Note this is done separately for "sine" and "cosine" components.
  
  nj = maxorder/2 + 1;
  for (m=0; m<mmax; m++) {
    uvec jsin(nj);
    uvec jcos(nj);
    mat Q(nrow, nj);
    mat R(nj, nj);
    vec sgn(nj);
    for (i=0; i<nj; i++) {
      jcos(i) = (m+i+1)*(m+i+1) - 2*m -1;
      jsin(i) = jcos(i) + 1;
    }
    
    qr_econ(Q, R, zm.cols(jcos));
    sgn = sign(R.diag());
    
    annzm.cols(jcos) = zpnorm * Q * diagmat(sgn);
    if (m > 0) {
      qr_econ(Q, R, zm.cols(jsin));
      sgn = sign(R.diag());
      annzm.cols(jsin) = zpnorm * Q * diagmat(sgn);
    }
    --nj;
  }
  
  //  highest order only has one term
  
  j = mmax*mmax;
  
  annzm.col(j) = zm.col(j) * zpnorm / norm(zm.col(j));
  annzm.col(j+1) = zm.col(j+1) * zpnorm / norm(zm.col(j+1));
  
  return annzm;
}
