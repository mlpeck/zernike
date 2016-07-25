// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rzernike
NumericVector rzernike(NumericVector rho, int n, int m);
RcppExport SEXP zernike_rzernike(SEXP rhoSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    __result = Rcpp::wrap(rzernike(rho, n, m));
    return __result;
END_RCPP
}
// zpmC
NumericMatrix zpmC(NumericVector rho, NumericVector theta, int maxorder);
RcppExport SEXP zernike_zpmC(SEXP rhoSEXP, SEXP thetaSEXP, SEXP maxorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type maxorder(maxorderSEXP);
    __result = Rcpp::wrap(zpmC(rho, theta, maxorder));
    return __result;
END_RCPP
}
