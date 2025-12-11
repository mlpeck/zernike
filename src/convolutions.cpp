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

//' 2D separable convolution - real matrix
//'
//' 2D convolution with a separable kernel
//'
//' Having a separable kernel implies that a
//' direct 2D convolution can be executed as a series
//' of 1D column- and row-wise convolutions.
//' This *may* be faster a direct 2D convlution or
//' using FFT's, and certainly benefits from multithreading.
//' The input matrix is implicitly zero padded, which differs
//' from an fft based convolution, where the input is implicitly
//' extended to infinity in all directions, therefore the outputs
//' will in general differ
//'
//' @param X a real valued matrix.
//' @param kernel
//'
//' @return a real matrix with the same dimension as X.
//[[Rcpp::export]]
mat conv2_sep_real(const mat& X, const vec& kernel){
  return conv2_sep(X, kernel);
}

//' 2D separable convolution - complex matrix
//'
//' 2D convolution with a separable kernel
//'
//' Having a separable kernel implies that a
//' direct 2D convolution can be executed as a series
//' of 1D column- and row-wise convolutions.
//' This *may* be faster a direct 2D convlution or
//' using FFT's, and certainly benefits from multithreading.
//' The input matrix is implicitly zero padded, which differs
//' from an fft based convolution, where the input is implicitly
//' extended to infinity in all directions, therefore the outputs
//' will in general differ
//'
//' @param X a complex valued matrix.
//' @param kernel
//'
//' @return a complex valued  matrix with the same dimension as X.
//[[Rcpp::export]]
cx_mat conv2_sep_complex(const cx_mat& X, const vec& kernel) {
  uword ksize = kernel.n_elem;
  cx_vec ckernel(ksize);
  ckernel.set_real(kernel);
  return conv2_sep(X, ckernel);
}

//' 2D Gaussian blur
//'
//' 2D Gaussian blur of a real matrix
//'
//' @param X a matrix
//' @param sigma standard deviation in pixels of the Gaussian kernel
//'
//' @return A matrix with the same dimensions as the input
//'
//' Performs direct convolution using [conv2_sep()] with a Gaussian kernel
//' sized to span 5 standard deviations.
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

//' 2D Gaussian blur of a complex valued matrix
//'
//' 2D Gaussian blur of a complex matrix
//'
//' @param X a complex valued matrix
//' @param sigma standard deviation in pixels of the Gaussian kernel
//'
//' @return A complex valued matrix with the same dimensions as the input
//'
//' Performs direct convolution using [conv2_sep_complex()] with a Gaussian kernel
//' sized to span 5 standard deviations. Preliminary benchmarking indicates this
//' is no faster than performing separate convolutions of the real and imaginary
//' components and combining the results.
//[[Rcpp::export]]
cx_mat gblur_complex(const cx_mat& X, const double sigma) {
  uword ksize = 5. * sigma;
  if (ksize % 2 == 0) ++ksize;
  vec kernel(ksize);
  cx_vec ckernel(ksize);
  vec y(ksize);
  y = arma::linspace(-2.5, 2.5, ksize);
  kernel = arma::normpdf(y);
  kernel /= arma::sum(kernel);
  ckernel.set_real(kernel);
  return conv2_sep(X, ckernel);
}

//' Derivative of Gaussian filter
//'
//' 2D filter that combines smoothing and (central) differencing.
//'
//' @param X a matrix
//' @param sigma standard deviation in pixels of the Gaussian kernel
//'
//' @return A named list with components Dx, Dy, Dxy, and the calculated kernel
//'
//' Calls arma::conv2(). Called by [circle.auto()].
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

//' 2D convolution - real matrix
//'
//' General 2D dircet convolution
//'
//' General 2D direct convolution with an arbitrary matrix kernel.
//' This is a wrapper for the Armadillo function `conv1`.
//'
//' @param X a real valued matrix.
//' @param kernel a real valued matrix
//'
//' @return a real matrix with the same dimension as X.
//[[Rcpp::export]]
mat convolve2d(const mat& X, const mat& kernel) {
  return arma::conv2(X, kernel, "same");
}
