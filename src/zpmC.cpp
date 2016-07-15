//
//	Copyright (C) 2016 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: MIT <https://opensource.org/licenses/MIT>
//

// Fill a matrix with Zernike polynomial values

# include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix zpmC(NumericVector rho, NumericVector theta, int maxorder) {

  int i, m, n, nn, n0, np2, nm2, order, nmax, mmax= maxorder/2;
  int nr = rho.size();
  int ncol = (mmax+1)*(mmax+1);
  double a0, a1, a2;
  double cosmtheta[mmax], sinmtheta[mmax];
  NumericMatrix zm(nr, ncol);
  
  //do some rudimentary error checking
  
  if (rho.size() != theta.size()) stop("Numeric vectors must be same length");
  if ((maxorder % 2) != 0) stop("maxorder must be even");
  
  //good enough

  
  for (n=0; n<nr; n++) {
      
      //cache values of cos and sin
      for (m=0; m<mmax; m++) {
          cosmtheta[m] = cos((m+1)*theta[n]);
          sinmtheta[m] = sin((m+1)*theta[n]);
      }
      
      zm(n, 0) = 1.0;                     //piston term
      zm(n, 1) = rho[n];                  //tilt
      zm(n, 3) = 2.*rho[n]*rho[n] - 1.;   //defocus
      
      // now fill in columns with m=n for n>0
      
      for (m=2; m <= mmax; m++) {
          zm(n, m*m) = pow(rho[n], (double) m);
      }
      
      //now the rest of the radially symmetric terms
      
      for (i=2; i<= mmax; i++) {
          nn = 2*i-2;
          np2 = (i*i+2*i);
          nm2 = (i*i-2*i);
          n0 =  (i*i-1);
          a0 = 1./(nn+2.);
          a1 = 4.*(nn+1);
          a2 = 2.*(nn+1);
          zm(n, np2) = a0*((a1*rho[n]*rho[n]-a2)*zm(n, n0)-nn*zm(n, nm2));
      }
      
      // non-symmetric terms
      
      for (m=1; m<mmax; m++) {
          order = 2*m+2;
          n0 = m*m;
          np2 = n0 + order +1;
          a0 = (double) (m+2.)/(4.*(m+1.)); //coefficients simplify for R(n=m+2,m)
          a1 = (double) 4.*(m+1.);          //in terms of R(n=m, m)
          a2 = (double) 4.*m + 4./(m+2.);
          zm(n, np2) = a0*(a1*rho[n]*rho[n] - a2)*zm(n, n0); //R(n=m+2,m) in terms of R(n=m,m)
          
          nmax = maxorder-m;
          if (nmax < (m+4)) break;
          for (nn = m+2; nn < nmax; nn += 2) {
              order += 2;
              nm2 = n0;
              n0 = np2;
              np2 = n0 + order + 1;
              a0 = (double) (nn+2.)/((nn+2.)*(nn+2.)-m*m);
              a1 = 4.*(nn+1);
              a2 = (double) (nn+m)*(nn+m)/nn+(nn-m+2.)*(nn-m+2.)/(nn+2.);
              zm(n, np2) = a0*((a1*rho[n]*rho[n] - a2)*zm(n, n0)
              - (double) (nn*nn-m*m)/nn * zm(n, nm2));
          }
      }
      
      // now multiply each column by normalizing factor and cos, sin
      
      n0 = 1;
      for (order=2; order <= maxorder; order+=2) {
          for(m=order/2; m>0; m--) {
              nn=order-m;
              a0 = sqrt(2.*(nn+1));
              zm(n, n0+1) = a0*sinmtheta[m-1]*zm(n, n0);
              zm(n, n0) *= a0*cosmtheta[m-1];
              n0 += 2;
          }
          nn = order;
          zm(n, n0) *= sqrt(nn+1.);
          n0++;
      }
  }
  return zm;
}
