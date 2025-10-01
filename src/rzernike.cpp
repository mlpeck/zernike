//
//	Copyright (C) 2016 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// radial Zernike polynomials
// Iterative version of the recurrence relation

// [[Rcpp::depends(RcppArmadillo)]]

# include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]

vec rzernike(const vec& rho, const int& n, const int& m) {
    
    if (((n-m) % 2) != 0 || n<0 || m<0 || n<m) Rcpp::stop("Bad index");

    uword np=rho.n_elem;
    vec rpoly(np);
    
    
    if ((n==0) && (m==0)) {
        rpoly.ones();
        return rpoly;
    } 
    
    if (n == m) {
        return pow(rho, m);
    } 
    
    if ((n==2) && (m==0)) {
        return 2.*rho % rho - 1;
    }
    
    vec r2(np), rp2(np), rm2(np);
    double a, a2, a4;
        
    r2 = rho % rho;
    rpoly = pow(rho, m);
    rm2.zeros();
    for (int j=m; j<=n-2; j+=2) {
        a = (double) ((j+2.)*(j+2.) - m*m)/(j+2.);
        a2 = (double) (j-m+2.)*(j-m+2.)/(j+2.);
        a4 = 0.;
        if (j > 0) {
            a2 += (double) (j+m)*(j+m)/j;
            a4 = (double) (m*m - j*j)/j;
        }
        
        rp2 = ((4.*(j+1.)*r2 - a2) % rpoly + a4 * rm2)/a;
        rm2 = rpoly;
        rpoly = rp2;
    }
    return rpoly;
}

        
