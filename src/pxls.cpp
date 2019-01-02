// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]

mat pxls(const mat& im, const vec& delta, const mat& tilt,
                    const vec& x, const vec& y) {
  uword nr = im.n_rows;
  uword nf = im.n_cols;

  colvec ph(nf);
  mat A(nf, 3);
  colvec z(nf);
  colvec b(3);
  mat B(nr, 3);
  
  for (uword n=0; n<nr; n++) {
    ph = delta + 4.0 * M_PI * (tilt.col(0) * x(n) + tilt.col(1) * y(n));
    A = join_rows(join_rows(ones(nf), cos(ph)), sin(ph));
    z = im.row(n).t();
    b = pinv(A) * z;
    B.row(n) = b.t();
  }
  return B;
}

  
