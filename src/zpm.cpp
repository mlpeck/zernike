//
//	Copyright (C) 2019 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// Fill a matrix with Zernike polynomial values

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <omp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]

mat zpm(const vec& rho, const vec& theta, const int& maxorder = 14) {

  int mmax = maxorder/2;
  int ncol = (mmax+1)*(mmax+1);
  uword imax = rho.size();
  mat zm(imax, ncol);
  
  
  //do some rudimentary error checking
  
  if (rho.size() != theta.size()) stop("Numeric vectors must be same length");
  if ((maxorder % 2) != 0) stop("maxorder must be even");
  
  //good enough
  
  #pragma omp parallel for
  for (uword i=0; i<imax; i++) {
    
    vec cosmtheta(mmax), sinmtheta(mmax);
    
    
    //cache values of cos and sin
    cosmtheta(0) = std::cos(theta(i));
    sinmtheta(0) = std::sin(theta(i));
    for (int m=1; m<mmax; m++) {
      cosmtheta(m) = cosmtheta(m-1)*cosmtheta(0) - sinmtheta(m-1)*sinmtheta(0);
      sinmtheta(m) = sinmtheta(m-1)*cosmtheta(0) + cosmtheta(m-1)*sinmtheta(0);
    }
    
    zm(i, 0) = 1.0;                     //piston term
    zm(i, 3) = 2. * rho(i) * rho(i) - 1.; //defocus
    
    // now fill in columns with m=n for n>0
    
    for (int m=1; m <= mmax; m++) {
      zm(i, m*m) = rho(i) * zm(i, (m-1)*(m-1));
    }

    // non-symmetric terms
    
    int order, n, m, nm, nm1mm1, nm1mp1, nm2m;
    
    for (order=4; order<=maxorder; order+=2) {
      for (m=order/2-1; m>0; m--) {
        n=order-m;
        nm = order*order/4 + n - m;
        nm1mm1 = (order-2)*(order-2)/4 + n - m;
        nm1mp1 = nm - 2;
        nm2m = nm1mm1 - 2;
        zm(i, nm) = rho(i)*(zm(i, nm1mm1) + zm(i, nm1mp1)) - zm(i, nm2m);
      }
      
      // m=0 (symmetric) term
      nm = order*order/4 + order;
      nm1mp1 = nm-2;
      nm2m = (order-2)*(order-2)/4+order-2;
      zm(i, nm) = 2.*rho(i)*zm(i, nm1mp1) - zm(i, nm2m);
    }

    // now multiply each column by normalizing factor and cos, sin
    
    int n0 = 1;
    double a0;
    for (order=2; order <= maxorder; order+=2) {
      for(m=order/2; m>0; m--) {
        n=order-m;
        a0 = sqrt(2.*(n+1));
        zm(i, n0+1) = a0*sinmtheta(m-1)*zm(i, n0);
        zm(i, n0) *= a0*cosmtheta(m-1);
        n0 += 2;
      }
      n = order;
      zm(i, n0) *= sqrt(n+1.);
      n0++;
    }
  }
  return zm;
}
