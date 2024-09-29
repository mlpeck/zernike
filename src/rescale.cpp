// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

//[[Rcpp::export]]

mat rescale(const mat& img, const double scale) {
  uword nr = img.n_rows;
  uword nc = img.n_cols;
  
  uword nr_s, nc_s;
  
  nr_s = scale * nr;
  nc_s = scale * nc;
  
  mat img_s(nr_s, nc_s);
  
  vec x_in = linspace(0., 1., nc);
  vec y_in = linspace(0., 1., nr);
  vec x_s = linspace(0., 1., nc_s);
  vec y_s = linspace(0., 1., nr_s);
  
  interp2(x_in, y_in, img, x_s, y_s, img_s);
  return(img_s);
}
