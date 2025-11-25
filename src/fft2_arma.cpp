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

// [[Rcpp::export]]
cx_mat fft_pad(const mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

// [[Rcpp::export]]
cx_mat fft_cx_pad(const cx_mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

// [[Rcpp::export]]
cx_mat fft(const mat& X) {
  return arma::fft2(X);
}

// [[Rcpp::export]]
cx_mat fft_cx(const cx_mat& X) {
  return arma::fft2(X);
}

// [[Rcpp::export]]
cx_mat ifft(const cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  return arma::ifft2(X);
}


// [[Rcpp::export]]
mat ifft_real(const cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  mat Y(nr, nc);
  Y = real(arma::ifft2(X));
  return Y;
}
