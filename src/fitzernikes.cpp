// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <omp.h>

using namespace arma;
using namespace Rcpp;

mat zapm(const vec&, const vec&, const double&, const int&, const int&);
mat zapm_iso(const vec&, const vec&, const double&, const int&, const int&);
mat zapm_128(const vec&, const vec&, const double&, const int&, const int&);
mat zapm_iso_128(const vec&, const vec&, const double&, const int&, const int&);
mat zpm(const vec&, const vec&, const int&);
mat zpm_cart(const vec&, const vec&, const int&, const bool&);

// [[Rcpp::export]]

vec fitzernikes(const vec& wf, const vec& rho, const vec& theta,
                const double& eps = 0.0, int maxorder = 14,
                int nthreads = -1, bool isoseq = false,
                bool usecirc = false, bool ext_prec = false) {
  
  uword nr = wf.size();
  uword ncol;
  
  if (rho.size() != theta.size()) Rcpp::stop("Coordinate vectors must be same length");
  if (rho.size() != nr) Rcpp::stop("Coordinate vectors must be same length as wf");
  
  if (nthreads <= 0) {
    nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);
  }
  
  if (isoseq) {
    ncol = (maxorder + 1) * (maxorder + 2) /2;
  } else {
    ncol = (maxorder/2 + 1) * (maxorder/2 + 1);
  }
  mat zm(nr, ncol);
  
  if (!isoseq) {
    if ((maxorder % 2) != 0) Rcpp::stop("maxorder must be even");
    
    if (eps <= 0.0 | usecirc) {
      zm = zpm(rho, theta, maxorder);
    } else {
      int nq = 2*maxorder + 5;
      if (ext_prec) {
        zm = zapm_128(rho, theta, eps, maxorder, nq);
      } else {
        zm = zapm(rho, theta, eps, maxorder, nq);
      }
    }
  } else {
   
    if (eps <= 0.0 | usecirc) {
      zm = zpm_cart(rho % cos(theta), rho % sin(theta), maxorder, true);
    } else {
      int nq = 2*maxorder + 5;
      if (ext_prec) {
        zm = zapm_iso_128(rho, theta, eps, maxorder, nq);
      } else {
        zm = zapm_iso(rho, theta, eps, maxorder, nq);
      }
    }
  }
  vec fit(nr);;
  
  fit = solve(zm.t() * zm, zm.t() * wf, solve_opts::fast + solve_opts::likely_sympd);
  
  return fit;
}
