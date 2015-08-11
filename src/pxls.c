#include <math.h>
#include <lapacke.h>

void pxls(int *nf, int *nr, double *im, double *delta,
		  double *tx, double *ty, double *x, double *y, double *B, int *info) {
	int M = *nf, N = *nr, m, n, i;
	lapack_int lda = 3, ldb = 1, nrhs=1, ncol = 3;
	double *im0 = im, *b0 = B;
	double A[M][3], work[M], ph;
	double *a0 = &A[0][0];
	
	for (n=0; n<N; n++) {
		for (m=0; m<M; m++) {
			ph = delta[m] + 4.*M_PI*(tx[m]*x[n]+ty[m]*y[n]);
			A[m][0] = 1.;
			A[m][1] = cos(ph);
			A[m][2] = sin(ph);
			work[m] = *im0++;
		}
		*info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N', (lapack_int) M, ncol,
		  nrhs, a0, lda, work, ldb);
		for (i=0; i<3; i++) {
			*b0 = work[i];
			++b0;
		}
	}
}

