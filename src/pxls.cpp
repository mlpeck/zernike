// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

//[[Rcpp::export]]

arma::mat pxls(arma::mat im, arma::vec delta, arma::mat tilt,
                    arma::vec x, arma::vec y) {
  int nr = im.n_rows;
  int nf = im.n_cols;
  int i, m, n;
  double ph;
  arma::mat B(nr, 3);
  
  arma::mat A(nf, 3);
  arma::colvec z(nf);
  arma::colvec b(3);
  
  for (n=0; n<nr; n++) {
    for (m=0; m<nf; m++) {
      ph = delta(m) + 4.0*M_PI*(tilt(m, 0)*x(n)+tilt(m, 1)*y(n));
      A(m, 0) = 1.;
      A(m, 1) = cos(ph);
      A(m, 2) = sin(ph);
      z(m) = im(n, m);
    }
    b = arma::solve(A, z);
    for (i=0; i<3; i++) {
      B(n, i) = b(i);
    }
  }
  return B;
}

  
