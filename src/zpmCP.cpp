//
//	Copyright (C) 2016-2019 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// Fill a matrix with Zernike polynomial values

// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <mutex>

using namespace Rcpp;
using namespace RcppParallel;

inline void zpm_row(const double rho, const double theta, const int maxorder,
                    RMatrix<double>::Row zm_row) {
  std::size_t m, n, n0, mmax = maxorder/2;
  std::size_t order, nm, nm1mm1, nm1mp1, nm2m;
  double a0;
  std::mutex mu;
  std::vector<double> cosmtheta(mmax);
  std::vector<double> sinmtheta(mmax);
  std::lock_guard<std::mutex> l(mu);

  //cache values of cos and sin
  cosmtheta[0] = std::cos(theta);
  sinmtheta[0] = std::sin(theta);
  for (m=1; m<mmax; m++) {
    cosmtheta[m] = cosmtheta[m-1]*cosmtheta[0] - sinmtheta[m-1]*sinmtheta[0];
    sinmtheta[m] = sinmtheta[m-1]*cosmtheta[0] + cosmtheta[m-1]*sinmtheta[0];
  }
  
  zm_row[0] = 1.0;                     //piston term
  zm_row[3] = 2.*rho*rho - 1.;   //defocus
  
    // now fill in columns with m=n for n>0
    
    for (m=1; m <= mmax; m++) {
      zm_row[m*m] = rho * zm_row[(m-1)*(m-1)];
    }

    // non-symmetric terms
    
    for (order=4; order<=maxorder; order+=2) {
      for (m=order/2-1; m>0; m--) {
        n=order-m;
        nm = order*order/4 + n - m;
        nm1mm1 = (order-2)*(order-2)/4 + n - m;
        nm1mp1 = nm - 2;
        nm2m = nm1mm1 - 2;
        zm_row[nm] = rho*(zm_row[nm1mm1] + zm_row[nm1mp1]) - zm_row[nm2m];
      }
      
      // m=0 (symmetric) term
      nm = order*order/4 + order;
      nm1mp1 = nm-2;
      nm2m = (order-2)*(order-2)/4+order-2;
      zm_row[nm] = 2.*rho*zm_row[nm1mp1] - zm_row[nm2m];
    }

  // now multiply each column by normalizing factor and cos, sin
  
  n0 = 1;
  for (order=2; order <= maxorder; order+=2) {
    for(m=order/2; m>0; m--) {
      n=order-m;
      a0 = std::sqrt(2.*(n+1));
      zm_row[n0+1] = a0*sinmtheta[m-1]*zm_row[n0];
      zm_row[n0] = a0*cosmtheta[m-1]*zm_row[n0];
      n0 += 2;
    }
    n = order;
    zm_row[n0] = zm_row[n0] * std::sqrt(n+1.);
    n0++;
  }
};  

struct FillZPM: public Worker{
  const RVector<double> rho;
  const RVector<double> theta;
  const std::size_t maxorder;
  RMatrix<double> zm;
  
  FillZPM(const NumericVector rho, const NumericVector theta, const int maxorder, 
          NumericMatrix zm)
    : rho(rho), theta(theta), maxorder(maxorder), zm(zm) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t n=begin; n<end; n++) {
      RMatrix<double>::Row zm_row = zm.row(n);
      zpm_row(rho[n], theta[n], maxorder, zm_row);
    }
  }
};

// [[Rcpp::export]]

NumericMatrix zpmCP(const NumericVector& rho, const NumericVector& theta, const int& maxorder = 14) {
  
  //do some rudimentary error checking
  
  if (rho.size() != theta.size()) stop("Numeric vectors must be same length");
  if ((maxorder % 2) != 0) stop("maxorder must be even");
  
  //good enough

  std::size_t mmax = maxorder/2;
  std::size_t ncol = (mmax+1) * (mmax+1);
  std::size_t nr = rho.size();
  
  NumericMatrix zm(nr, ncol);
  
  FillZPM fillzpm(rho, theta, maxorder, zm);
  parallelFor(0, rho.size(), fillzpm);
  return zm;
}
