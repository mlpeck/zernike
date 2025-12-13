// [[Rcpp::depends(RcppArmadillo)]]

/*******************
 *
 * 2d fftshift and ifftshift
 *
 * *****************/

# include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' 2D fftshift
//'
//' Shift quadrants of a matrix
//'
//' This swaps the quadrants of
//' a matrix to put the (0,0) element
//' in the center. The most common use is
//' to center the DC component of a Fourier Transformed matrix
//' for visualization or to simplify other matrix operations
//' @param X a real or complex valued matrix.
//' @return a complex matrix with the same dimension as X.
//' @seealso [ifftshift()] is the inverse operation.
//' @examples
//' X <- matrix(1:16, 4, 4)
//' XS <- fftshift(X)
//' XS
//' ifftshift(XS)
// [[Rcpp::export]]
cx_mat fftshift(cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  cx_mat X_buf(nr, nc);

  uword nr2 = nr/2;
  uword nc2 = nc/2;

  if (nr % 2) {
    X_buf.head_rows(nr2)   = X.tail_rows(nr2);
    X_buf.tail_rows(nr2+1) = X.head_rows(nr2+1);
  } else {
    X_buf.head_rows(nr2) = X.tail_rows(nr2);
    X_buf.tail_rows(nr2) = X.head_rows(nr2);
  }

  if (nc % 2) {
    X.head_cols(nc2)   = X_buf.tail_cols(nc2);
    X.tail_cols(nc2+1) = X_buf.head_cols(nc2+1);
  } else {
    X.head_cols(nc2) = X_buf.tail_cols(nc2);
    X.tail_cols(nc2) = X_buf.head_cols(nc2);
  }

  return X;
}

//' 2D ifftshift
//'
//' Inverse fftShift
//'
//' Inverts the fftshift operation,
//' swapping quadrants of the input matrix.
//'
//' @param X a real or complex valued matrix.
//' @return a complex matrix with the same dimension as X.
//' @seealso [fftshift()].
//' @examples
//' X <- matrix(1:16, 4, 4)
//' XS <- fftshift(X)
//' XS
//' ifftshift(XS)
// [[Rcpp::export]]
cx_mat ifftshift(cx_mat& X) {
  uword nr = X.n_rows;
  uword nc = X.n_cols;
  cx_mat X_buf(nr, nc);

  uword nr2 = nr/2;
  uword nc2 = nc/2;

  if (nr % 2) {
    X_buf.head_rows(nr2+1) = X.tail_rows(nr2+1);
    X_buf.tail_rows(nr2)   = X.head_rows(nr2);
  } else {
    X_buf.head_rows(nr2) = X.tail_rows(nr2);
    X_buf.tail_rows(nr2) = X.head_rows(nr2);
  }

  if (nc % 2) {
    X.head_cols(nc2+1) = X_buf.tail_cols(nc2+1);
    X.tail_cols(nc2)   = X_buf.head_cols(nc2);
  } else {
    X.head_cols(nc2) = X_buf.tail_cols(nc2);
    X.tail_cols(nc2) = X_buf.head_cols(nc2);
  }

  return X;
}
