/**********************
 *
 * 2d FFT's wrapping FFTW
 *
 * FFTW is Copyright © 2003, 2007-11 Matteo Frigo, Copyright © 2003, 2007-11 Massachusetts Institute of Technology.
 *
 * This portion is licensed under GPL, version 2 or later per terms of FFTW license.
 *
 * ref: Matteo Frigo and Steven G. Johnson, “The design and implementation of FFTW3,” Proc. IEEE 93 (2), 216–231 (2005).
 *
 * ********************/


// [[Rcpp::depends(RcppArmadillo)]]

# include <RcppArmadillo.h>

using namespace arma;

#include <vector>
#include <complex>
#include <cmath>
#include <armadillo>
#include <fftw3.h>


static fftw_plan plan_fw {nullptr};
static fftw_plan plan_bw {nullptr};
static uword nr_last {0};
static uword nc_last {0};


/************
 *
 * 2D fft with real input matrix, complex out
 *
 * **********/

//' 2D FFT - FFTW3 based version
//'
//' 2D FFT wrapping routines in fftw3
//'
//' @param X a real valued matrix
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//' This is mostly an Armadillo based C++ wrapper to functions in fftw3 On first
//' call or if matrix dimensions change calls `fftw_plan_dft_2d` with the planning flag
//' `FFTW_MEASURE` (currently hard coded in the source), then `fftw_execute_dft` to perform the transform.
//' This gives a reasonable balance between overhead of the first call and performance
//' improvement in subsequent calls compared to the Armadillo based fft routine.
//' Nevertheless this is probably only useful if dozens to hundreds (or more) frames are
//' to be processed.
//'
//' @seealso [fft()], [fft_fftw_cx()], [ifft_fftw()]
//'
//' @references Matteo Frigo and Steven G. Johnson, "The Design and Implementation of FFTW3,"
//' roceedings of the IEEE 93 (2), 216–231 (2005). Invited paper, Special Issue on Program Generation, Optimization, and Platform Adaptation.
// [[Rcpp::export]]
cx_mat fft_fftw(const mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;

  cx_mat XC(nr, nc);
  XC.set_real(X);
  cx_vec xcvec = XC.as_col();
  cx_vec planvec = xcvec;

  if (plan_fw == nullptr || nr != nr_last || nc != nc_last) {
    plan_fw = fftw_plan_dft_2d(
      nc, nr,
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      FFTW_FORWARD, FFTW_MEASURE);
    nr_last = nr;
    nc_last = nc;
  }

  fftw_execute_dft(plan_fw,
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()),
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()));

  XC = reshape(xcvec, nr, nc);
  return XC;
}

/************
 *
 * 2D fft with complex input matrix, complex out
 *
 * **********/

//' 2D FFT - complex valued input -FFTW3 based version
//'
//' 2D FFT wrapping routines in fftw3
//'
//' @param X a complex valued matrix
//'
//' @returns A complex valued matrix with the discrete Fourier transformed input.
//'
//'
//' @seealso [fft_cx()], [fft_fftw()], [ifft_fftw()]
// [[Rcpp::export]]
cx_mat fft_fftw_cx(const cx_mat& XC) {
  uword nr = XC.n_rows;
  uword nc = XC.n_cols;

  cx_vec xcvec = XC.as_col();
  cx_mat XCT(nr, nc);
  cx_vec planvec = xcvec;

  if (plan_fw == nullptr || nr != nr_last || nc != nc_last) {
    plan_fw = fftw_plan_dft_2d(
      nc, nr,
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      FFTW_FORWARD, FFTW_MEASURE);
    nr_last = nr;
    nc_last = nc;
  }

  fftw_execute_dft(plan_fw,
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()),
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()));

  XCT = reshape(xcvec, nr, nc);
  return XCT;
}

//' 2D Inverse FFT - FFTW3 based version
//'
//' 2D Inverse FFT wrapping routines in fftw3
//'
//' @param X a complex valued matrix
//'
//' @returns A complex valued matrix with the scaled inverse Fourier transformed input.
//'
//'
//' @seealso [ifft()], fft_fftw()]
// [[Rcpp::export]]
cx_mat ifft_fftw(const cx_mat& XC) {
  uword nr = XC.n_rows;
  uword nc = XC.n_cols;

  cx_vec xcvec = XC.as_col();
  cx_vec planvec = xcvec;
  cx_mat XCT(nr, nc);

  if (plan_bw == nullptr|| nr != nr_last || nc != nc_last) {
    plan_bw = fftw_plan_dft_2d(
      nc, nr,
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      reinterpret_cast<fftw_complex*>(planvec.memptr()),
      FFTW_BACKWARD, FFTW_MEASURE);
    nr_last = nr;
    nc_last = nc;
  }

  fftw_execute_dft(plan_bw,
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()),
                   reinterpret_cast<fftw_complex*>(xcvec.memptr()));

  XCT = reshape(xcvec, nr, nc)/(nr*nc);
  return XCT;
}

