
// [[Rcpp::depends(RcppArmadillo)]]

#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>

using namespace Rcpp;
using namespace arma;

mat pwrap(const mat& phase) {
  uword nr = phase.n_rows;
  uword nc = phase.n_cols;
  uword ne = phase.n_elem;
  mat wphase(nr, nc);
  
  for (uword i=0; i<ne; i++) {
    wphase(i) = fmod(phase(i) + M_PI, 2.*M_PI) - M_PI;
  }
  
  return wphase;
}

/*********
 * 
 * PSI by least squares with optional per frame weights
 * 
*/

//[[Rcpp::export]]

List lspsiC(const mat& images, const rowvec& phases, const vec& wt) {
  int nf = images.n_cols;
  int np = images.n_rows;
  
  if (phases.n_elem != nf || wt.n_elem != nf) stop("Dimension mismatch");
  mat S = join_cols(join_cols(ones<rowvec>(nf), cos(phases)), sin(phases));
  
  mat B = images * diagmat(sqrt(wt)) * pinv(S * diagmat(sqrt(wt)));
  vec phi = atan2(-B.col(2), B.col(1));
  vec mod = sqrt(square(B.col(1)) + square(B.col(2)));
  mod = mod/max(mod);
  
  return List::create(Named("phi") = phi,
                      Named("mod") = mod);
}

/***************
 * The "Advanced iterative algorithm of Wang & Han 2004
 * As interpreted by ML Peck
 * License: MIT
*/


//[[Rcpp::export]]

List aiapsiC(const mat& images, const rowvec& phases_init, const double& ptol, const int& maxiter, const bool& trace) {
  uword M = images.n_rows;
  uword N = images.n_cols;
  
  rowvec phases = phases_init - phases_init(0);
  rowvec phases_last = phases_init - phases_init(0);
  
  mat S(3, N);
  mat Phi(M, 3);
  vec phi(M);
  vec mod(M);
  double sdp;
  int i;
  
  for (i=0; i<maxiter; i++) {

    // Estimate the phase from the current estimate of phase shifts
    
    S = join_cols(join_cols(ones<rowvec>(N), cos(phases)), sin(phases));
    Phi = images * pinv(S);
    phi = atan2(-Phi.col(2), Phi.col(1));
    Phi = join_rows(join_rows(ones(M), cos(phi)), -sin(phi));
    
    // Estimate phase shifts from current estimate of phase.
    // Solves normal equations for speed
    
    S = pinv(Phi.t() * Phi) * Phi.t() * images;
    phases = atan2(S.row(2), S.row(1));
    phases = phases - phases(0);
    sdp = norm(sin(phases-phases_last), 2);
    
    if (trace) {
      cout << "iteration: "<< i << "sdp = " << sdp << " Phases = " << phases;
    }
    
    // repeat until convergence
    
    if (sdp < ptol) break;
    phases_last = phases;
  }
  
  // one more calculation of the phase
  
  S = join_cols(join_cols(ones<rowvec>(N), cos(phases)), sin(phases));
  Phi = images * pinv(S);
  phi = atan2(-Phi.col(2), Phi.col(1));
  mod = sqrt(square(Phi.col(1)) + square(Phi.col(2)));
  mod = mod/max(mod);
  return List::create(Named("phi") = phi,
                      Named("mod") = mod,
                      Named("phases") = phases,
                      Named("iter") = i);
}

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
  mat Phi(M, 3);
  mat V(N, 3);
  vec phases(N), phases_last(N);
  mat P(3, 3), P_last(3, 3);
  mat S(N, 3);
  mat vp(N, 3);
  double sdp, normp;
  
  svd_econ(U_full, d_full, V_full, images);
  V = V_full.cols(0, 2);
  Phi = U_full.cols(0, 2) * diagmat(d_full.subvec(0, 2));
  phases_last = atan2(V.col(2), V.col(1));
  phases_last =  phases_last - phases_last(0);
  S = join_rows(join_rows(ones(N), cos(phases_last)), sin(phases_last));
  P = V.t() * S;
  P_last = P;
  
  for (int i=0; i<maxiter; i++) {
    vp = V * P;
    phases = atan2(vp.col(2), vp.col(1));
    phases = phases - phases(0);
    S = join_rows(join_rows(ones(N), cos(phases)), sin(phases));
    P = V.t() * S;
    sdp = norm(sin(phases-phases_last), 2);
    normp = norm(P - P_last, 2);
    if (trace) {
      cout << "Iteration: " << i << " sdp = " << sdp << " normp = " << normp << " phases = " << phases;
    }
    if (sdp < ptol || normp < ptol) break;
    phases_last = phases;
    P_last = P;
  }
  Phi = Phi * pinv(P).t();
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



/********************
 * 
 * the c++ part of tiltpsi
 * 
*/


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

  
