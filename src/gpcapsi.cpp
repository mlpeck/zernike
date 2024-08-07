
// [[Rcpp::depends(RcppArmadillo)]]

#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>

using namespace Rcpp;
using namespace arma;


/*************
 * 
 * My generalized PCA algorithm
 * 
*/

//[[Rcpp::export]]

List gpcapsiC(const mat& images, const double& ptol, const int& maxiter, const bool& trace) {
  uword M = images.n_rows;
  uword N = images.n_cols;
    
  mat U_full(M, N);
  mat V_full(N, N);
  vec d_full(N);
  mat U_cov(N, N);
  mat Phi(M, 3);
  mat V(N, 3);
  rowvec phases(N), phases_last(N);
  mat P(3, 3), P_last(3, 3);
  mat S(3, N);
  mat pv(3, N);
  double sdp, normp;
  
  // estimate phase shifts from pixelwise mean subtracted image data
  
  mat img_sub(M, N);
  colvec img_mean(M);
  mat DV(N, 2);

  img_sub = images;
  img_mean = mean(images, 1);  
  img_sub.each_col() -= img_mean; 
  
  svd_econ(U_full, d_full, V_full, img_sub);
  DV = V_full.cols(0, 1) * diagmat(d_full.subvec(0, 1));
  phases_last = atan2(-DV.col(1), DV.col(0)).as_row();
  phases_last =  phases_last - phases_last(0);
  cout << "Initial phase shift guess" << phases_last << endl;
  S = join_cols(ones<rowvec>(N), cos(phases_last), sin(phases_last));

  svd_econ(U_full, d_full, V_full, images);
  V = V_full.cols(0, 2);
  Phi = U_full.cols(0, 2) * diagmat(d_full.subvec(0, 2));
  
  P = S * V;
  P_last = P;
  
  for (int i=0; i<maxiter; i++) {
    pv = P * V.t();
    phases = atan2(pv.row(2), pv.row(1));
    phases = phases - phases(0);
    S = join_cols(ones<rowvec>(N), cos(phases), sin(phases));
    P = S * V;
    sdp = norm(sin(phases-phases_last), 2);
    normp = norm(P - P_last, 2);
    if (trace) {
      cout << "Iteration: " << i << " sdp = " << sdp << " normp = " << normp << " phases = " << phases;
    }
    if (sdp < ptol || normp < ptol) break;
    phases_last = phases;
    P_last = P;
  }
  Phi = Phi * pinv(P);
  vec phi = atan2(-Phi.col(2), Phi.col(1));
  vec mod = sqrt(square(Phi.col(1)) + square(Phi.col(2)));
  mod = mod/max(mod);
  double r2 = sum(square(d_full.subvec(0, 2)))/sum(square(d_full));
  r2 = sqrt(r2/(1.-r2));
  return List::create(Named("phi") = phi,
                      Named("mod") = mod,
                      Named("phases") = phases,
                      Named("P") = P,
                      Named("snr") = r2,
                      Named("eigen") = square(d_full));
}

