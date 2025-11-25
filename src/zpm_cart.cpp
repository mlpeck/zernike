/**********************
 *
 * Copyright © 2022 Michael Peck <mlpeck54 -at- gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 **************************/



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]


/*****************
 * 
 * Convert a matrix of Zernike polynomial values
 * from unit scaled to unit variance (orthonormal scaling).
 * Columns must be in ISO/ANSI sequence. 
 * Only check performed is on number of columns.
 * 
*****************/


//' Normalize matrix of Zernike polynomial values.
//'
//' Convert a matrix of Zernike polynomial values from
//' unit scaled to unit variance aka orthonormal form.
//'
//' @param uzpm matrix of Zernike polynomial values
//' @param maxorder the maximum radial order.
//'
//' @return matrix in orthonormal form.
//'
//' @details
//'  This is intended only for ISO/ANSI ordered matrices. The
//'  only check performed is that the number of columns in the
//'  matrix matches the expected number given by the argument
//'  `maxorder`.
//'  This is called by [gradzpm()] and [zpm_cart()]
//'  if `unit_variance` is set to `true` in the respective
//'  function calls.
//' @md
// [[Rcpp::export]]
mat norm_zpm(mat& uzpm, const int& maxorder) {
  uword ncol = uzpm.n_cols;
  int n, m, mp, j=0;
  double zmult;
  
  if (ncol != (maxorder+1)*(maxorder+2)/2) Rcpp::stop("Matrix is wrong size for entered order");
  
  for (n=0; n<=maxorder; n++) {
    for (m=0; m<=n; m++) {
      mp = n - 2*m;
      zmult = std::sqrt(n+1.0);
      if (mp != 0) {
        zmult *= sqrt(2.0);
      }
      uzpm.col(j) = zmult * uzpm.col(j);
      ++j;
    }
  }
  return uzpm;
}
      
  
  
/********************
 Fill matrixes with:
 unit scaled Zernikes in cartesian coordinates
 x and y directional derivatives
 
 Note: This returns ZP values in ISO?ANSI ordering with sine terms on
       the "left" side of the pyramid.
 
 Source:
 Anderson, T.B. (2018) Optics Express 26, #5, 18878
 URL: https://doi.org/10.1364/OE.26.018878 (open access)
********************/


//' Zernike polynomials and cartesian gradients
//'
//' Calculate Zernike polynomial values and Cartesian gradients in
//' ISO/ANSI sequence for a set of Cartesian coordinates.
//'
//' @param x a vector of x coordinates for points on a unit disk.
//' @param y a vector of y coordinates.
//' @param maxorder the maximum radial polynomial order (defaults to 12).
//' @param unit_variance logical: return with orthonormal scaling? (default `false`)
//' @param return_zpm logical: return Zernike polynomial matrix? (default `true`)
//'
//' @return a named list with the matrices `zm` (optional but returned by default), `dzdx`, `dzdy`.
//'
//' @references
//'   Anderson, T.B. (2018) Optics Express 26, #5, 18878
//'   <https://doi.org/10.1364/OE.26.018878> (open access)
//'
//' @details
//'  Uses the recurrence relations in the above publication to calculate Zernike
//'  polynomial values and their directional derivatives in Cartesian coordinates. These are
//'  known to be both efficient and numerically stable.
//'
//'  Columns are in ISO/ANSI sequence: for each radial order n >= 0 the azimuthal orders m are sequenced
//'  m = {-n, -(n-2), ..., (n-2), n}, with sine components for negative m and cosine for positive m. Note this
//'  is the opposite ordering from the extended Fringe set and the ordering of aberrations is quite different. 
//'  For example the two components of trefoil are in the 7th and 10th column while coma is in
//'  columns 8 and 9 (or 7 and 8 with 0-indexing). Note also that except for tilt and coma-like aberrations
//'  (m=1) non-axisymmetric aberrations will be separated.
//'
//'  All three matrices will have the same dimensions on return. Columns 0 and 1 of `dzdx` will be all 0,
//'  while columns 0 and 3 of `dzdy` are 0.
//'
//' @seealso [zpm()] uses the same recurrence relations for polar coordinates and extended
//'  Fringe set ordering, which is the more common indexing scheme for optical design/testing
//'  software.
//' @seealso [zpm_cart()] calculates and returns the Zernike polynomial values only.
//'
//' @examples
//'  rho <- seq(0.2, 1., length=5)
//'  theta <- seq(0, 1.6*pi, length=5)
//'  rt <- expand.grid(theta, rho)
//'  x <- c(0, rt[,2]*cos(rt[,1]))
//'  y <- c(0, rt[,2]*sin(rt[,1]))
//'  gzpm <- gradzpm(x, y)
//'
//' @md
// [[Rcpp::export]]
List gradzpm(const vec& x, const vec& y, const int& maxorder = 12, 
             const bool& unit_variance = true, const bool& return_zpm = true) {
  
  int n, m;
  int ncol = (maxorder+1)*(maxorder+2)/2;
  uword nrow = x.size();
  int j;
  bool n_even;
  mat zm(nrow, ncol), dzdx(nrow, ncol), dzdy(nrow, ncol);
  
  //do some rudimentary error checking
  
  if (x.n_elem != y.n_elem) stop("Numeric vectors must be same length");
  if (maxorder < 1) stop("maxorder must be >= 1");
  
  // good enough
  
  
  // starting values for recursions
  
  zm.col(0).ones();
  zm.col(1) = y;
  zm.col(2) = x;
  
  dzdx.col(0).zeros();
  dzdx.col(1).zeros();
  dzdx.col(2).ones();
  
  dzdy.col(0).zeros();
  dzdy.col(1).ones();
  dzdy.col(2).zeros();
  
  j = 3;
  n_even = true;
  for (n=2; n <= maxorder; n++) {
    
    // fill in  m=0
    m=0;
    zm.col(j) = x % zm.col(j-n) + y % zm.col(j-1);
    dzdx.col(j) = (double) n * zm.col(j-n);
    dzdy.col(j) = (double) n * zm.col(j-1);
    ++j;
    
    for (m=1; m < n/2; m++) {
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      dzdx.col(j) = (double) n*zm.col(j-n) + (double) n*zm.col(j-n-1) + dzdx.col(j-2*n);
      dzdy.col(j) = (double) n*zm.col(j-2*m-1) - (double) n*zm.col(j-2*m) + dzdy.col(j-2*n);
      ++j;
    }
    if (n_even) {
      zm.col(j) = 2.*x % zm.col(j-n) + 2.*y % zm.col(j-n-1) - zm.col(j-2*n);
      dzdx.col(j) = 2.*(double) n*zm.col(j-n) + dzdx.col(j-2*n);
      dzdy.col(j) = 2.*(double) n*zm.col(j-2*m-1) + dzdy.col(j-2*n);
      ++m;
      ++j;
    } else {
      zm.col(j) = y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      dzdx.col(j) = (double) n*zm.col(j-n-1) + dzdx.col(j-2*n);
      dzdy.col(j) = (double) n*(zm.col(j-2*m-1) - zm.col(j-2*m)) + dzdy.col(j-2*n);
      ++m;
      ++j;
      
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - zm.col(j-2*n);
      dzdx.col(j) = (double) n*(zm.col(j-n) + zm.col(j-n-1)) + dzdx.col(j-2*n);
      dzdy.col(j) = (double) n*zm.col(j-2*m-1) + dzdy.col(j-2*n);
      ++m;
      ++j;
    }
    for (m=m; m<n; m++) {
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      dzdx.col(j) = (double) n*zm.col(j-n) + (double) n*zm.col(j-n-1) + dzdx.col(j-2*n);
      dzdy.col(j) = (double) n*zm.col(j-2*m-1) - (double) n*zm.col(j-2*m) + dzdy.col(j-2*n);
      ++j;
    }
    
    // m = n
    
    zm.col(j) = x % zm.col(j-n-1) - y % zm.col(j-2*n);
    dzdx.col(j) = (double) n * zm.col(j-n-1);
    dzdy.col(j) = (double) -n * zm.col(j-2*n);
    ++j;
    
    n_even = !n_even;
  }
  
  //orthonormalize if requested
    
  if (unit_variance) {
    zm = norm_zpm(zm, maxorder);
    dzdx = norm_zpm(dzdx, maxorder);
    dzdy = norm_zpm(dzdy, maxorder);
  }
  
  if (return_zpm) {
    return List::create(Named("zm") = zm, Named("dzdx") = dzdx, Named("dzdy") = dzdy);
  }
  return List::create(Named("dzdx") = dzdx, Named("dzdy") = dzdy);
}

/*****************
 * Fill matrix of Zernike polynomial values in cartesian coordinates.
 * 
 * This is just the above code with all references to directional
 * derivatives removed.
*****************/


//' Zernike polynomials
//'
//' Calculate Zernike polynomial values in
//' ISO/ANSI sequence for a set of Cartesian coordinates.
//'
//' @param x a vector of x coordinates for points on a unit disk.
//' @param y a vector of y coordinates.
//' @param maxorder the maximum radial polynomial order (defaults to 12).
//' @param unit_variance logical: return with orthonormal scaling? (default `false`)
//'
//' @return a matrix of Zernike polynomial values evaluated at the input
//'  Cartesian coordinates and all radial and azimuthal orders from
//'  0 through `maxorder`.
//'
//' @references
//'   Anderson, T.B. (2018) Optics Express 26, #5, 18878
//'   <https://doi.org/10.1364/OE.26.018878> (open access)
//'
//' @details This is the same algorithm and essentially the same code as [gradzpm()]
//'  except directional derivatives aren't calculated.
//' @md
// [[Rcpp::export]]
mat zpm_cart(const vec& x, const vec& y, const int& maxorder = 14, 
                       const bool& unit_variance = true) {
  
  int n, m;
  int ncol = (maxorder+1)*(maxorder+2)/2;
  uword nrow = x.n_elem;
  int j;
  bool n_even;
  mat zm(nrow, ncol);
  
  //do some rudimentary error checking
  
  if (x.n_elem != y.n_elem) stop("Numeric vectors must be same length");
  if (maxorder < 1) stop("maxorder must be >= 1");
  
  // good enough
  
  
  // starting values for recursions
  
  zm.col(0).ones();
  zm.col(1) = y;
  zm.col(2) = x;
  
  
  j = 3;
  n_even = true;
  for (n=2; n <= maxorder; n++) {
    
    // fill in  m=0
    m=0;
    zm.col(j) = x % zm.col(j-n) + y % zm.col(j-1);
    ++j;
    
    for (m=1; m < n/2; m++) {
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      ++j;
    }
    if (n_even) {
      zm.col(j) = 2.*x % zm.col(j-n) + 2.*y % zm.col(j-n-1) - zm.col(j-2*n);
      ++m;
      ++j;
    } else {
      zm.col(j) = y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      ++m;
      ++j;
      
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - zm.col(j-2*n);
      ++m;
      ++j;
    }
    for (m=m; m<n; m++) {
      zm.col(j) = x % zm.col(j-n) + y % zm.col(j-2*m-1) + x % zm.col(j-n-1) - y % zm.col(j-2*m) - zm.col(j-2*n);
      ++j;
    }
    
    // m = n
    
    zm.col(j) = x % zm.col(j-n-1) - y % zm.col(j-2*n);
    ++j;
    
    n_even = !n_even;
  }
  
  //orthonormalize if requested
    
  if (unit_variance) {
    zm = norm_zpm(zm, maxorder);
  }
  
  return zm;
}

// [[Rcpp::export]]
vec lsfit_qr(const vec& wf, const mat& zpm) {
  return solve(zpm, wf);
}

// [[Rcpp::export]]
vec lsfit_norm(const vec& wf, const mat& zpm) {
  return solve(zpm.t() * zpm, zpm.t() * wf);
}

