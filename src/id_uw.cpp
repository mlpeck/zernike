#include <Rcpp.h>
#include "puw.h"

using namespace Rcpp;

#define unwrap(dp) {          \
    puw[kn] = puw[kh] + dp;   \
    flags[kn] |= UNWRAPPED;   \
    uw[kn] = MASK;            \
    todo[ntodo++] = kn;       \
}

//' Compiled code via Rcpp for Itoh's method of phase unwrapping
//' 
//' Called by [brcutpuw()] for fast phase unwrapping
//' 
//' @param nr number of rows in phase matrix
//' @param nc number of columns in phase matrix
//' @param phase phase matrix converted to vector
//' @param mask matrix of mask values converted to vector
//' @param dx wrapped phase differences in x direction
//' @param dy wrapped phase differences in y direction
//' @return a vector with the unwrapped phase
//' @details
//'   This is called by [brcutpuw()] through [idiffpuw()] 
//'   but is also user callable.
//'   Wrapped phase values and differences are divided by \code{2*pi} before input
//'   making the input values in the range [-1/2, 1/2).
//'   In [brcutpuw()] the mask indicates areas outside the interferogram area
//'   and lines of branch cuts
//' @seealso [brcutpuw()], [idiffpuw()]
//' @author M.L. Peck (mlpeck54 -at- gmail.com)  


// [[Rcpp::export]]

NumericVector id_dxy_uw(const int& nr, const int& nc, 
                        const NumericVector& phase, const NumericVector& mask,
                        const NumericVector& dx, const NumericVector& dy,
                        IntegerVector uw) {
  size_t len = phase.length();
  NumericVector puw(len);
  std::vector<char> flags(len);
  std::vector<size_t> todo(len);
  size_t kh, kn, x, y;
  size_t ntodo = 0;
  
  for (size_t k=0; k<len; k++) flags[k] = traits::is_na<REALSXP>(mask[k]);
  for (size_t k=0; k<len; k++) {
    if (!flags[k]) {
      kn = kh = k;
      break;
    }
  }
  puw[kh] = phase[kh];
  unwrap(0.0);
  
  while (ntodo) {
    kh = todo[--ntodo];
    x = kh % nr;
    y = kh / nr;
    kn = kh-1;
    if (x > 0 && !flags[kn])
      unwrap(-dx[kn]);
    kn = kh+1;
    if (x < (nr-1) && !flags[kn])
      unwrap(dx[kh]);
    kn = kh - nr;
    if (y > 0 && !flags[kn])
      unwrap(-dy[kn]);
    kn = kh + nr;
    if (y < (nc-1) && !flags[kn])
      unwrap(dy[kh]);
  }
  return puw;
}

// same routine as above, but without dx and dy passed.

#define unwrap2(dp) {         \
    puw[kn] = puw[kh] + dp;   \
    flags[kn] |= UNWRAPPED;   \
    todo[ntodo++] = kn;       \
}

//' Compiled code via Rcpp for Itoh's method of phase unwrapping
//' 
//' Called by [idiffpuw()] for fast phase unwrapping
//' 
//' @param nr number of rows in phase matrix
//' @param nc number of columns in phase matrix
//' @param phase phase matrix converted to vector
//' @return a vector with the unwrapped phase
//' @details
//'   This is called by [idiffpuw()] 
//'   but is also user callable.
//'   Wrapped phase values are divided by \code{2*pi} before input
//'   making the input values in the range [-1/2, 1/2).
//'   In [brcutpuw()] the mask indicates areas outside the interferogram area
//'   and lines of branch cuts
//' @seealso [brcutpuw()], [idiffpuw()]
//' @author M.L. Peck (mlpeck54 -at- gmail.com)  


// [[Rcpp::export]]

NumericVector id_uw(const int& nr, const int& nc, const NumericVector& phase) {
  size_t len = phase.length();
  NumericVector puw(len);
  std::vector<char> flags(len);
  std::vector<size_t> todo(len);
  size_t kh, kn, x, y;
  size_t ntodo = 0;
  
  for (size_t k=0; k<len; k++) flags[k] = traits::is_na<REALSXP>(phase[k]);
  for (size_t k=0; k<len; k++) {
    if (!flags[k]) {
      kn = kh = k;
      break;
    }
  }
  puw[kh] = phase[kh];
  unwrap2(0.0);
  
  while (ntodo) {
    kh = todo[--ntodo];
    x = kh % nr;
    y = kh / nr;
    kn = kh-1;
    if (x > 0 && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh+1;
    if (x < (nr-1) && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh - nr;
    if (y > 0 && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
    kn = kh + nr;
    if (y < (nc-1) && !flags[kn])
      unwrap2(WRAP(phase[kn]-phase[kh]));
  }
  return puw;
}

