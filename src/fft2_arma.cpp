/**********************
 *
 * Copyright © 2026 Michael Peck <mlpeck54 -at- gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 **************************/


// [[Rcpp::depends(RcppArmadillo)]]

/*******************
 *
 * 2d FFT's using only Armadillo supplied routines
 *
 * *****************/

# include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' 2D FFT
//'
//' 2D FFT wrapping Armadillo's `fft2`
//'
//' @param X a real valued matrix
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//' This is just a wrapper for Armadillo's `fft2`
//' function. It works well for dimensions that are
//' *highly composite*, not just powers of 2 as
//' implied by the online documentation. `X` need not
//' be a square matrix.
//' @references Conrad Sanderson and Ryan Curtin.
//' Armadillo: An Efficient Framework for Numerical Linear Algebra.
//' International Conference on Computer and Automation Engineering, 2025.
//' <br>
//' Conrad Sanderson and Ryan Curtin.
//' Practical Sparse Matrices in C++ with Hybrid Storage and Template-Based Expression Optimisation.
//' Mathematical and Computational Applications, Vol.24, No.3, 2019.
// [[Rcpp::export]]
cx_mat fft(const mat& X) {
  return arma::fft2(X);
}

//' 2D FFT - complex valued input
//'
//' 2D FFT wrapping Armadillo's `fft2`
//'
//' @param X a complex valued matrix
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//' @seealso [fft()]
// [[Rcpp::export]]
cx_mat fft_cx(const cx_mat& X) {
  return arma::fft2(X);
}

//' 2D FFT with padding
//'
//' 2D FFT wrapping Armadillo's `fft2`
//'
//' @param X a real valued matrix
//' @param npad size to pad to
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//' This is provided for cases where the input matrix has one or both dimensions
//' that are *not* highly composite (that is, have prime factors that are larger than 5).
//' `npad` itself should be highly composite. This isn't checked. The R function [nextn()]
//' can be used to find a suitable size.
//'
//' @seealso [fft()]
// [[Rcpp::export]]
cx_mat fft_pad(const mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

//' 2D FFT with padding - complex valued input
//'
//' 2D FFT wrapping Armadillo's `fft2`
//'
//' @param X a complex valued matrix
//' @param npad size to pad to
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//' This is provided for cases where the input matrix has one or both dimensions
//' that are *not* highly composite (that is, have prime factors that are larger than 5).
//' `npad` itself should be highly composite. This isn't checked. The R function [nextn()]
//' can be used to find a suitable size.
//'
//' @seealso [fft()]
// [[Rcpp::export]]
cx_mat fft_cx_pad(const cx_mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

//' 2D Inverse FFT - complex valued input
//'
//' 2D Inverse FFT wrapping Armadillo's `ifft2`
//'
//' @param X a complex valued matrix
//'
//' @returns A complex valued matrix with the inverse discrete Fourier transformed input.
//'
//' In contrast to some implementations this returns the *scaled* inverse transform.
//'
//' @seealso [fft()]
// [[Rcpp::export]]
cx_mat ifft(const cx_mat& X) {
  return arma::ifft2(X);
}

//' 2D Inverse FFT - real valued output
//'
//' 2D Inverse FFT wrapping Armadillo's `ifft2`
//'
//' @param X a complex valued matrix
//'
//' @returns A real valued matrix with the real component of the inverse discrete Fourier transformed input.
//'
//' In many cases where a real valued matrix is transformed, some operations performed in the Fourier domain,
//' and then inverse transformed the output will be real but with a small imaginary component.
//' This function simply discards the imaginary component. No check is made that this is
//' actually appropriate.
//'
//' @seealso [fft()]
// [[Rcpp::export]]
mat ifft_real(const cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  mat Y(nr, nc);
  Y = real(arma::ifft2(X));
  return Y;
}

//' 2D Gaussian blur
//'
//' 2D Gaussian blur -- FFT based version
//'
//' @param X a real valued matrix
//' @param sigma standard deviation in pixels of the Gaussian kernel
//'
//' @returns the smoothed input matrix
//'
//' Numerical experiments suggest that it's almost always faster to use
//' the direct convolution based Gaussian blur in [gblur()].
// [[Rcpp::export]]
mat gblur_fft(const mat& X, const double sigma) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  uword ksize = 5. * sigma;
  if (ksize % 2 == 0) ++ksize;
  vec kernel(ksize);
  vec y(ksize);
  y = arma::linspace(-2.5, 2.5, ksize);
  kernel = arma::normpdf(y);
  kernel /= arma::sum(kernel);
  mat kmat = kernel * kernel.t();
  return ifft_real(arma::fft2(X) % arma::fft2(kmat, nr, nc));
}

