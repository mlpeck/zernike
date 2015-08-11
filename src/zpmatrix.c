//
//	Copyright (C) 2012 Michael Peck <mpeck1 -at- ix.netcom.com>
//
//	License: GPL <http://www.gnu.org/licenses/licenses.html#GPL>
//

// Fill a matrix with Zernike polynomial values

# include <math.h>

void zpmatrix(double *rho, double *theta, int *nr, int *maxorder, double *zm) {
  int i, m, n, nn, order, nmax, mmax= *maxorder/2;
  int c0, cp2, cm2;
  long s0, sm2, sp2;
  double a0, a1, a2;

  // fill in column 0 (n=0, m=0)

  for (n=0; n<(*nr); n++) zm[n] = 1.0;

  // now fill in columns with m=n for n>0

  for (m=1; m <= mmax; m++) {
	s0 = m*m*(*nr);
	for (n=0; n<(*nr); n++) zm[s0++] = pow(rho[n], (double) m);
  }

  //defocus

  s0 = 3*(*nr);
  for (n=0; n<(*nr); n++) zm[s0++] = 2*rho[n]*rho[n] -1.0;

  //now the rest of the radially symmetric terms

  for (i=2; i<= mmax; i++) {
	nn= 2*i-2;
	sp2=(i*i+2*i)*(*nr);
	s0 = (i*i-1)*(*nr);
	sm2 = (i*i-2*i)*(*nr);
	a0 = 1./(nn+2.);
	a1 = 4.*(nn+1);
	a2 = 2.*nn+2.;
	for (n=0; n<(*nr); n++) 
	  zm[sp2++] = a0*((a1*rho[n]*rho[n]-a2)*zm[s0++]-nn*zm[sm2++]);
  }

  // get the rest

  for (m=1; m<mmax; m++) {
	nn = m;
	order = 2*m+2;
	c0 = m*m;
	s0 = c0*(*nr);
	cp2 = c0 + order +1;
	sp2 = cp2*(*nr);
	a0 = (double) (nn+2.)/((nn+2.)*(nn+2.)-m*m);
	a1 = (double) 4.*(nn+1.);
	a2 = (double) (nn+m)*(nn+m)/nn+(nn-m+2.)*(nn-m+2.)/(nn+2.);
	for (n=0; n<(*nr); n++) {
	  zm[sp2++] = 
	   a0*(a1*rho[n]*rho[n] - a2)*zm[s0++];
	}
	nmax = *maxorder-m;
	if (nmax < (m+4)) break;
	for (nn = m+2; nn < nmax; nn+=2) {
	  order += 2;
	  cm2 = c0;
	  c0 = cp2;
	  cp2 = c0 + order + 1;
	  sm2 = cm2*(*nr);
	  s0 = c0 *(*nr);
	  sp2 = cp2*(*nr);
	  a0 = (double) (nn+2.)/((nn+2.)*(nn+2.)-m*m);
	  a1 = (double) 4.*(nn+1.);
	  a2 = (double) (nn+m)*(nn+m)/nn+(nn-m+2.)*(nn-m+2.)/(nn+2.);
	  for (n=0; n<(*nr); n++) {
		zm[sp2++] = 
		a0*((a1*rho[n]*rho[n] -a2)*zm[s0++]
		-(double) (nn*nn-m*m)/nn * zm[sm2++]);
	  }
	}
  }

  // now multiply each column by normalizing factor and cos, sin

  c0 = 1;
  for (order=2; order <= (*maxorder); order+=2) {
	for(m=order/2; m>0; m--) {
	  nn=order-m;
	  s0 = c0*(*nr);
	  a0 = sqrt(2.*(nn+1));
	  for (n=0; n<(*nr); n++) {
		zm[s0+(*nr)] = a0*sin(m*theta[n])*zm[s0];
		zm[s0] = a0*cos(m*theta[n])*zm[s0];
		s0++;
	  }
	  c0 += 2;
	}
	nn = order;
	s0 = c0*(*nr);
	a0 = sqrt(nn+1.);
	for (n=0; n<(*nr); n++) zm[s0++] = a0*zm[s0];
	c0++;
  }
}

  