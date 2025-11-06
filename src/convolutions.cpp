// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]

using namespace arma;

template<typename Tmat, typename Tvec>
Tmat conv2_sep(const Tmat& X, const Tvec& kernel) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  uword nk = kernel.size();

  Tmat X_int(nr, nc), X_out(nr, nc);
  #pragma omp parallel for
  for (uword i = 0; i < nc; ++i) {
    X_int.col(i) = arma::conv(X.col(i), kernel, "same");
  }

  #pragma omp parallel for
  for (uword i = 0; i < nr; ++i) {
    X_out.row(i) = arma::conv(X_int.row(i), kernel.t(), "same");
  }

  return X_out;
}


//[[Rcpp::export]]
mat conv2_sep_real(const mat& X, const vec& kernel){
  return conv2_sep(X, kernel);
}

//[[Rcpp::export]]
cx_mat conv2_sep_complex(const cx_mat& X, const vec& kernel) {
  uword ksize = kernel.n_elem;
  cx_vec ckernel(ksize);
  ckernel.set_real(kernel);
  return conv2_sep(X, ckernel);
}

//[[Rcpp::export]]
mat gblur(const mat& X, const double sigma) {
  uword ksize = 5. * sigma;
  if (ksize % 2 == 0) ++ksize;
  vec kernel(ksize);
  vec y(ksize);
  y = arma::linspace(-2.5, 2.5, ksize);
  kernel = arma::normpdf(y);
  kernel /= arma::sum(kernel);
  return conv2_sep(X, kernel);
}


//[[Rcpp::export]]
Rcpp::List d_of_g(const mat& X, const double sigma) {
  uword nr, nc;
  nr = X.n_rows;
  nc = X.n_cols;
  mat Dx(nr, nc), Dy(nr, nc), Dxy(nr, nc);

  uword ksize = 5. * sigma;
  if (ksize % 2 == 0) ++ksize;
  vec x(ksize), xval(ksize);
  rowvec y(ksize), yval(ksize);
  mat kernel(ksize, ksize);

  x = arma::linspace(-2.5, 2.5, ksize);
  xval = x % arma::normpdf(x);

  y = arma::linspace<rowvec>(-2.5, 2.5, ksize);
  yval = arma::normpdf(y);

  kernel = xval * yval;
  Dx = arma::conv2(X, kernel, "same");
  Dy = arma::conv2(X, kernel.t(), "same");
  Dxy = arma::square(Dx) + arma::square(Dy);

  return Rcpp::List::create(Rcpp::Named("Dx") = Dx,
                            Rcpp::Named("Dy") = Dy,
                            Rcpp::Named("Dxy") = Dxy,
                            Rcpp::Named("kernel") = kernel);
}


//[[Rcpp::export]]
mat convolve2d(const mat& X, const mat& kernel) {
  return arma::conv2(X, kernel, "same");
}
