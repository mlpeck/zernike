
// [[Rcpp::depends(RcppArmadillo)]]

#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>

using namespace Rcpp;
using namespace arma;

typedef struct {vec x; bool convergence; int iter;} lmreturn;

/**************
 * 
 * A stupidly simple implementation of a Levenberg-Marquardt
 * algorithm for nonlinear least squares.
 * 
 * Adapted directly from the algorithm description on p. 27
 * of the booklet "Methods for Non-Linear Least Squares
 * Problems" by Madsen, Nielsen & Tingleff (2004)
 * download from http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 * 
*/

lmreturn mylevmar(vec (*res_fun)(const vec&, const List&),
             mat (*jac_fun)(const vec&, const List&),
             uword m, int n, vec x0, List adata,
             int itmax = 100, vec control= {1.e-3, 2., sqrt(datum::eps), sqrt(datum::eps)}) {
  
  int k=0;
  double tau = control(0);
  double nu = control(1);
  double eps1 = control(2);
  double eps2 = control(3);
  
  mat A(n, n);
  vec g(n);
  vec feval(m);
  vec fnew(m);
  mat Jac(m, n);
  mat In = eye<mat>(n, n);
  vec x = x0;
  vec hx(n), xnew(n);
  bool found;
  double mu;
  double F, Fnew, dL, rho;
  lmreturn retvals;
  
  Jac = jac_fun(x, adata);
  feval = res_fun(x, adata);
  A = Jac.t() * Jac;
  g = Jac.t() * feval;
  mu = tau * max(A.diag());
  found = (norm(g, "Inf") <= eps1);
  
  while (!found && k<itmax) {
    k++;
    hx = -solve(A + mu*In, g);
    if (norm(hx, 2) <= eps2 * (norm(x, 2) + eps2)) {
      found = true;
      break;
    }
    xnew = x + hx;
    F = 0.5 * dot(feval, feval);
    fnew = res_fun(xnew, adata);
    Fnew = 0.5 * dot(fnew, fnew);
    dL = 0.5 * dot(hx, (mu*hx - g));
    rho = (F-Fnew)/dL;
    if (rho > 0.) {
      x = xnew;
      Jac = jac_fun(x, adata);
      A = Jac.t() * Jac;
      g = Jac.t() * fnew;
      feval = fnew;
      mu *= fmax(1./3., 1 - pow(2.*rho -1., 3));
      nu = control(1);
      found = (norm(g, "Inf") <= eps1);
    } else {
      mu *= nu;
      nu *= 2.;
    }
  }
  
  retvals.x = x;
  retvals.convergence = found;
  retvals.iter = k;
  
  return retvals;
  
}

// pixel-wise least squares with known phases, tilts, etc.

mat pxls(const mat& im, const rowvec& phases, const mat& zcs, const mat& coords) {
  uword nr = im.n_rows;
  uword nf = im.n_cols;

  rowvec ph(nf);
  mat A(3, nf);
  rowvec b(3);
  mat B(nr, 3);
  
  for (uword n=0; n<nr; n++) {
    ph = phases + coords.row(n) * zcs;
    A = join_cols(join_cols(ones<rowvec>(nf), cos(ph)), sin(ph));
    b = im.row(n) * pinv(A);
    B.row(n) = b;
  }
  return B;
}

// residuals and analytic Jacobian for mylevmar. Exported just in case I want to use minpack.lm

// [[Rcpp::export]]

vec res_frame(const vec& pars, const List& adata) {
  vec img = adata["img"];
  mat coords = adata["coords"];
  double abar = adata["abar"];
  double bbar = adata["bbar"];
  vec phi = adata["phi"];
  int np = pars.n_elem;
  vec ph(img.n_elem);
  ph = pars[0] + coords * pars.tail(np-1);
  return img - abar - bbar * cos(phi + ph);
}

// [[Rcpp::export]]

mat jac_frame(const vec& pars, const List& adata) {
  vec img = adata["img"];
  mat coords = adata["coords"];
  double bbar = adata["bbar"];
  vec phi = adata["phi"];
  int np = pars.n_elem;
  uword m = img.n_elem;
  vec ph(m);
  vec sp(m);
  mat jac(m, np);
  ph = phi + pars[0] + coords * pars.tail(np-1);
  sp = sin(ph);
  jac = join_rows(ones(m), coords);
  jac.each_col() %= sp;
  return bbar*jac;
}

/************
 * 
 * wrap phase into [-pi, pi)
 * This is slower than R version
 * but export it anyway with a different name.
 * 
*/

//[[Rcpp::export]]

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



// The main routine. This is what gets called from R

//[[Rcpp::export]]

List tiltpsiC(const mat& images, const rowvec& phases_init, const mat& coords,
              const double& ptol, const int& maxiter, const bool& trace) {
  int N = images.n_cols;
  uword M = images.n_rows;
  int np = coords.n_cols + 1;
  
  if (np < 3) stop("It doesn't make sense to use this algorithm without at least variable tilts.\nUse AIA or PC instead");
  if (phases_init.n_elem != N) stop("# phases must match # of images");
  if (coords.n_rows != M) stop("coordinates must have same # rows as image matrix");

  rowvec phases = phases_init - phases_init(0);
  
  mat S(3, N);
  mat Phi(M, 3);
  vec phi(M);
  vec mod(M);
  vec res(M);
  mat zcs = zeros<mat>(np-1, N);
  vec t0(np-1);
  vec sse = zeros<vec>(maxiter);
  double abar, bbar;
  List adata;
  adata["coords"] = coords;
  vec pars(np);
  mat pt_last = join_cols(phases, zcs);
  mat pt = pt_last;
  double dpt;
  int i;
  lmreturn lmrets;

  S = join_cols(join_cols(ones<rowvec>(N), cos(phases)), sin(phases));
  Phi = images * pinv(S);
  
  for (i=0; i<maxiter; i++) {
    abar = mean(Phi.col(0));
    bbar = mean(sqrt(square(Phi.col(1)) + square(Phi.col(2))));
    phi = atan2(-Phi.col(2), Phi.col(1));
    adata["abar"] = abar;
    adata["bbar"] = bbar;
    adata["phi"] = phi;
    
    for (int n=0; n<N; n++) {
      pars(0) = phases(n);
      pars.tail(np-1) = zcs.col(n);
      adata["img"] = images.col(n);
      lmrets = mylevmar((*res_frame), 
                      (*jac_frame), 
                      M, np, pars, adata);
      if (!lmrets.convergence) warning("lm convergence reported failed");
      pars = lmrets.x;
      res = res_frame(pars, adata);
      sse(i) = sse(i) + pow((double) norm(res, 2), 2);
      phases(n) = pars(0);
      zcs.col(n) = pars.tail(np-1);
    }
    phases = phases - phases(0);
    t0 = zcs.col(0);
    zcs.each_col() -= t0;
    pt = join_cols(phases, zcs);
    dpt = norm(pt - pt_last, 2);
    if (trace) {
      cout << "Iteration " << i << " sse = " << sse(i) << " dpt = " << dpt << endl;
    }
    if (dpt < ptol) break;
    Phi = pxls(images, phases, zcs, coords);
    pt_last = pt;
  }
  Phi = pxls(images, phases, zcs, coords);
  phi = atan2(-Phi.col(2), Phi.col(1));
  mod = sqrt(square(Phi.col(1)) + square(Phi.col(2)));
  mod = mod/max(mod);
  return List::create(Named("phi") = phi, 
                      Named("mod") = mod,
                      Named("phases") = pwrap(phases),
                      Named("zcs") = zcs/(2.*M_PI),
                      Named("iter") = i,
                      Named("sse") = sse);
}
