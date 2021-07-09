//
//	Copyright (C) 2016 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// radial Zernike polynomials
// Iterative version of the recurrence relation

# include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector rzernike(NumericVector rho, int n, int m) {
    
    unsigned int np=rho.size();
    unsigned int i;
    NumericVector rpoly(np);
    double a, a2, a4;
    
    if (((n-m) %2) != 0 || n<0 || m<0 || n<m) stop("Bad index");
    
    if ((n==0) && (m==0)) {
        rpoly.fill(1.0);
        return rpoly;
    } else if (n == m) {
        return pow(rho, m);
    } else if ((n==2) && (m==0)) {
        return 2.*rho*rho - 1;
    } else {
        NumericVector r2(np), rp2(np), rm2(np);
        int j;
        
        r2 = rho * rho;
        rpoly = pow(rho, m);
        rm2.fill(0.0);
        for (j=m; j<=n-2; j+=2) {
            a = (double) ((j+2.)*(j+2.) - m*m)/(j+2.);
            a2 = (double) (j-m+2.)*(j-m+2.)/(j+2.);
            a4 = 0.;
            if (j > 0) {
                a2 += (double) (j+m)*(j+m)/j;
                a4 = (double) (m*m - j*j)/j;
            }
            for (i=0; i<np; i++) {
                rp2[i] = ((4.*(j+1.)*r2[i] - a2) * rpoly[i] + a4 * rm2[i])/a;
                rm2[i] = rpoly[i];
                rpoly[i] = rp2[i];
            }
        }
        return rpoly;
    }
}

        
