#include <Rcpp.h>
#include "puw.h"

using namespace Rcpp;

#define SWAP(a, i, j)     \
{                         \
  size_t t = a[i];        \
  a[i] = a[j];            \
  a[j] = t;               \
}

void todo_push(const size_t ndx, const NumericVector& qual, std::vector<size_t>& todo, size_t& end) {
  size_t child;
  todo[end] = ndx;
  child = end++;
  while (child > 0) {
    size_t parent = (child-1)/2;
    if (qual[todo[parent]] < qual[todo[child]]) {
      SWAP(todo, parent, child);
      child = parent;
    } else {
      break;
    }
  }
}

size_t todo_pop(const NumericVector& qual, std::vector<size_t>& todo, size_t& end) {
  size_t result = todo[0];
  size_t root = 0;
  
  --end;
  SWAP(todo, 0, end);
  while(2*root+1 < end) {
    size_t child = 2*root+1;  //left child
    if (child+1 < end && qual[todo[child]] < qual[todo[child+1]]) {
      ++child;
    }
    if (qual[todo[root]] < qual[todo[child]]) {
      SWAP (todo, root, child);
      root = child;
    } else {
      break;
    }
  }
  return result;
}

// [[Rcpp::export]]

//' Compiled code via Rcpp for quality guided phase unwrapping
//' 
//' Called by [qpuw()] for fast quality guided phase unwrapping
//' 
//' @param nr number of rows in phase matrix
//' @param nc number of columns in phase matrix
//' @param phase phase matrix converted to vector
//' @param qual quality matrix converted to vector
//' @return a vector with the unwrapped phase
//' @details
//'   This is called by [qpuw()] but is also user callable.
//'   Wrapped phase values are divided by \code{2*pi} before input
//'   making the input values in the range [-1/2, 1/2).
//' @seealso [qpuw()], [idiffpuw()]
//' @author M.L. Peck (mlpeck54 -at- gmail.com)  
//'   with valuable programming advice from Steve Koehler
NumericVector q_uw(const int& nr, const int& nc, const NumericVector& phase, const NumericVector& qual) {
  size_t len = phase.length();
  NumericVector puw(len);
  std::vector<char> flags(len);
  std::vector<size_t> todo(len);
  
  double qm = -HUGE;
  size_t kh, kn, x, y, end=0;
  int nk[4] = {-1, 1, -nr, nr};
  int bxy[4] = {0, -(nr-1), 0, -(nc-1)};
  int cxy[4];
  
  for (size_t k=0; k<len; k++) {
    flags[k] = traits::is_na<REALSXP>(phase[k]);
    if (qual[k] > qm && !flags[k]) {
      kh = k;
      qm = qual[k];
    }
  };
  
  // start at highest quality point
  
  puw[kh] = phase[kh];
  flags[kh] |= UNWRAPPED;
  todo_push(kh, qual, todo, end);
  
  while (end) {
    kh = todo_pop(qual, todo, end);
    x = kh % nr;
    y = kh / nr;
    cxy[1] = -(cxy[0]=x);
    cxy[3] = -(cxy[2]=y);
    
    // neighbors of currently selected point
    
    for (int k=0; k<4; k++) {
      kn = kh+nk[k];
      if ((cxy[k] > bxy[k]) && !flags[kn]) {
        puw[kn] = puw[kh] + WRAP(phase[kn]-phase[kh]);
        flags[kn] |= UNWRAPPED;
        todo_push(kn, qual, todo, end);
      }
    }
  }
  return puw;
}
