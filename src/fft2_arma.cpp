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

//'  Wrappers for Armadillo 2D fft's with padding
//'
//'  2D ffts with real and complex inputs
//'
//'  @param X a real or complex valued metrix
//'  @param npad size to pad to (should be highly composite)
//'
//'  @returns A complex valued transformed matrix of size npad x npad.
//'
//'  @details For convenience there are versions for
//'    real or complex inputs. Armadillo's `fft2` function
//'    is poorly suited to matrixes of sizes other than
//'    highly composite. These functions zero pad the input
//'    and compute the transform of the padded matrix.
//'    `npad` obviously should be highly composite or a power of two.
// [[Rcpp::export]]
cx_mat fft_pad(const mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

//'  @rdname fft_pad
//'  @export
// [[Rcpp::export]]
cx_mat fft_cx_pad(const cx_mat& X, const uword npad) {
  if (npad < X.n_rows || npad < X.n_cols) {
    Rcpp::stop("npad must be larger than both matrix dimensions");
  }
  return arma::fft2(X, npad, npad);
}

//'  Wrappers for Armadillo 2D fft's
//'
//'  2D ffts with real and complex inputs
//'
//'  @param X a real or complex valued metrix
//'
//'  @returns A complex valued transformed matrix
//'
//'  @details For convenience there are versions for
//'    both real and complex inputs. These functions
//'    work well if the input dimensions are highly composite,
//'    that is have no factors larger than 5. If not,
//'    use the padded versions [fft_pad], [fft_cx_pad],
//'    or the FFTW3 based versions [fft_fftw], [fft_fftw_cx].
//'
//'  @references [Armadillo home page](https://arma.sourceforge.net/).
//'    citation: Conrad Sanderson and Ryan Curtin.
//'    *Armadillo: An Efficient Framework for Numerical Linear Algebra.
//'    International Conference on Computer and Automation Engineering*, 2025.
//'  @export
// [[Rcpp::export]]
cx_mat fft(const mat& X) {
  return arma::fft2(X);
}

//'  @rdname fft
//'  @export
// [[Rcpp::export]]
cx_mat fft_cx(const cx_mat& X) {
  return arma::fft2(X);
}

//'  Wrappers for Armadillo 2D inverse fft's
//'
//'  2D inverse ffts with complex inputs
//'
//'  @param X a complex valued metrix
//'
//'  @returns A complex or real valued transformed matrix
//'
//'  @details For convenience there are versions for
//'    both real and complex outputs. Note that a round trip
//'    real FFT -> inverse FFT will produce a numerically
//'    complex valued output with negligible imaginary component.
//'    The function [ifft_real()] is intended for this case.
//'    These functions
//'    work well if the input dimensions are highly composite,
//'    that is have no factors larger than 5.
//'  @export
// [[Rcpp::export]]
cx_mat ifft(const cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  return arma::ifft2(X);
}


//'  @rdname ifft
//'  @export
// [[Rcpp::export]]
mat ifft_real(const cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  mat Y(nr, nc);
  Y = real(arma::ifft2(X));
  return Y;
}
