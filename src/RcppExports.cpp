// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fitzernikes
vec fitzernikes(const vec& wf, const vec& rho, const vec& theta, const double& eps, int maxorder, int nthreads, bool isoseq, bool usecirc, bool ext_prec);
RcppExport SEXP _zernike_fitzernikes(SEXP wfSEXP, SEXP rhoSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP maxorderSEXP, SEXP nthreadsSEXP, SEXP isoseqSEXP, SEXP usecircSEXP, SEXP ext_precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type wf(wfSEXP);
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type isoseq(isoseqSEXP);
    Rcpp::traits::input_parameter< bool >::type usecirc(usecircSEXP);
    Rcpp::traits::input_parameter< bool >::type ext_prec(ext_precSEXP);
    rcpp_result_gen = Rcpp::wrap(fitzernikes(wf, rho, theta, eps, maxorder, nthreads, isoseq, usecirc, ext_prec));
    return rcpp_result_gen;
END_RCPP
}
// gpcapsiC
List gpcapsiC(const mat& images, const double& ptol, const int& maxiter, const bool& trace);
RcppExport SEXP _zernike_gpcapsiC(SEXP imagesSEXP, SEXP ptolSEXP, SEXP maxiterSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type images(imagesSEXP);
    Rcpp::traits::input_parameter< const double& >::type ptol(ptolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(gpcapsiC(images, ptol, maxiter, trace));
    return rcpp_result_gen;
END_RCPP
}
// id_dxy_uw
NumericVector id_dxy_uw(const int& nr, const int& nc, const NumericVector& phase, const NumericVector& mask, const NumericVector& dx, const NumericVector& dy, IntegerVector uw);
RcppExport SEXP _zernike_id_dxy_uw(SEXP nrSEXP, SEXP ncSEXP, SEXP phaseSEXP, SEXP maskSEXP, SEXP dxSEXP, SEXP dySEXP, SEXP uwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type phase(phaseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type dy(dySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type uw(uwSEXP);
    rcpp_result_gen = Rcpp::wrap(id_dxy_uw(nr, nc, phase, mask, dx, dy, uw));
    return rcpp_result_gen;
END_RCPP
}
// id_uw
NumericVector id_uw(const int& nr, const int& nc, const NumericVector& phase);
RcppExport SEXP _zernike_id_uw(SEXP nrSEXP, SEXP ncSEXP, SEXP phaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type phase(phaseSEXP);
    rcpp_result_gen = Rcpp::wrap(id_uw(nr, nc, phase));
    return rcpp_result_gen;
END_RCPP
}
// lspsiC
List lspsiC(const mat& images, const rowvec& phases, const vec& wt);
RcppExport SEXP _zernike_lspsiC(SEXP imagesSEXP, SEXP phasesSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type images(imagesSEXP);
    Rcpp::traits::input_parameter< const rowvec& >::type phases(phasesSEXP);
    Rcpp::traits::input_parameter< const vec& >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(lspsiC(images, phases, wt));
    return rcpp_result_gen;
END_RCPP
}
// aiapsiC
List aiapsiC(const mat& images, const rowvec& phases_init, const double& ptol, const int& maxiter, const bool& trace);
RcppExport SEXP _zernike_aiapsiC(SEXP imagesSEXP, SEXP phases_initSEXP, SEXP ptolSEXP, SEXP maxiterSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type images(imagesSEXP);
    Rcpp::traits::input_parameter< const rowvec& >::type phases_init(phases_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type ptol(ptolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(aiapsiC(images, phases_init, ptol, maxiter, trace));
    return rcpp_result_gen;
END_RCPP
}
// q_uw
NumericVector q_uw(const int& nr, const int& nc, const NumericVector& phase, const NumericVector& qual);
RcppExport SEXP _zernike_q_uw(SEXP nrSEXP, SEXP ncSEXP, SEXP phaseSEXP, SEXP qualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< const int& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type phase(phaseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type qual(qualSEXP);
    rcpp_result_gen = Rcpp::wrap(q_uw(nr, nc, phase, qual));
    return rcpp_result_gen;
END_RCPP
}
// readraw
NumericMatrix readraw(CharacterVector fname, NumericVector channels);
RcppExport SEXP _zernike_readraw(SEXP fnameSEXP, SEXP channelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type channels(channelsSEXP);
    rcpp_result_gen = Rcpp::wrap(readraw(fname, channels));
    return rcpp_result_gen;
END_RCPP
}
// rescale
mat rescale(const mat& img, const double scale);
RcppExport SEXP _zernike_rescale(SEXP imgSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type img(imgSEXP);
    Rcpp::traits::input_parameter< const double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rescale(img, scale));
    return rcpp_result_gen;
END_RCPP
}
// rzernike
vec rzernike(const vec& rho, const int& n, const int& m);
RcppExport SEXP _zernike_rzernike(SEXP rhoSEXP, SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rzernike(rho, n, m));
    return rcpp_result_gen;
END_RCPP
}
// pwrap
mat pwrap(const mat& phase);
RcppExport SEXP _zernike_pwrap(SEXP phaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type phase(phaseSEXP);
    rcpp_result_gen = Rcpp::wrap(pwrap(phase));
    return rcpp_result_gen;
END_RCPP
}
// tiltpsiC
List tiltpsiC(const mat& images, const rowvec& phases_init, const mat& coords, const double& ptol, const int& maxiter, const bool& trace);
RcppExport SEXP _zernike_tiltpsiC(SEXP imagesSEXP, SEXP phases_initSEXP, SEXP coordsSEXP, SEXP ptolSEXP, SEXP maxiterSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type images(imagesSEXP);
    Rcpp::traits::input_parameter< const rowvec& >::type phases_init(phases_initSEXP);
    Rcpp::traits::input_parameter< const mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const double& >::type ptol(ptolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(tiltpsiC(images, phases_init, coords, ptol, maxiter, trace));
    return rcpp_result_gen;
END_RCPP
}
// gol_welsch
vec gol_welsch(const double& eps, vec& qwts);
RcppExport SEXP _zernike_gol_welsch(SEXP epsSEXP, SEXP qwtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< vec& >::type qwts(qwtsSEXP);
    rcpp_result_gen = Rcpp::wrap(gol_welsch(eps, qwts));
    return rcpp_result_gen;
END_RCPP
}
// rzernike_ann
mat rzernike_ann(const vec& rho, const double& eps, const int& n, const int& m, const vec& xq, const vec& qwts);
RcppExport SEXP _zernike_rzernike_ann(SEXP rhoSEXP, SEXP epsSEXP, SEXP nSEXP, SEXP mSEXP, SEXP xqSEXP, SEXP qwtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const vec& >::type xq(xqSEXP);
    Rcpp::traits::input_parameter< const vec& >::type qwts(qwtsSEXP);
    rcpp_result_gen = Rcpp::wrap(rzernike_ann(rho, eps, n, m, xq, qwts));
    return rcpp_result_gen;
END_RCPP
}
// zapm
mat zapm(const vec& rho, const vec& theta, const double& eps, const int& maxorder, const int& nq);
RcppExport SEXP _zernike_zapm(SEXP rhoSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP maxorderSEXP, SEXP nqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const int& >::type nq(nqSEXP);
    rcpp_result_gen = Rcpp::wrap(zapm(rho, theta, eps, maxorder, nq));
    return rcpp_result_gen;
END_RCPP
}
// zapm_iso
mat zapm_iso(const vec& rho, const vec& theta, const double& eps, const int& maxorder, const int& nq);
RcppExport SEXP _zernike_zapm_iso(SEXP rhoSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP maxorderSEXP, SEXP nqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const int& >::type nq(nqSEXP);
    rcpp_result_gen = Rcpp::wrap(zapm_iso(rho, theta, eps, maxorder, nq));
    return rcpp_result_gen;
END_RCPP
}
// zapm_128
mat zapm_128(const vec& rho, const vec& theta, const double& eps, const int& maxorder, const int& nq);
RcppExport SEXP _zernike_zapm_128(SEXP rhoSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP maxorderSEXP, SEXP nqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const int& >::type nq(nqSEXP);
    rcpp_result_gen = Rcpp::wrap(zapm_128(rho, theta, eps, maxorder, nq));
    return rcpp_result_gen;
END_RCPP
}
// zapm_iso_128
mat zapm_iso_128(const vec& rho, const vec& theta, const double& eps, const int& maxorder, const int& nq);
RcppExport SEXP _zernike_zapm_iso_128(SEXP rhoSEXP, SEXP thetaSEXP, SEXP epsSEXP, SEXP maxorderSEXP, SEXP nqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const int& >::type nq(nqSEXP);
    rcpp_result_gen = Rcpp::wrap(zapm_iso_128(rho, theta, eps, maxorder, nq));
    return rcpp_result_gen;
END_RCPP
}
// zpm
mat zpm(const vec& rho, const vec& theta, const int& maxorder);
RcppExport SEXP _zernike_zpm(SEXP rhoSEXP, SEXP thetaSEXP, SEXP maxorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    rcpp_result_gen = Rcpp::wrap(zpm(rho, theta, maxorder));
    return rcpp_result_gen;
END_RCPP
}
// zpmCP
NumericMatrix zpmCP(const NumericVector& rho, const NumericVector& theta, const int& maxorder);
RcppExport SEXP _zernike_zpmCP(SEXP rhoSEXP, SEXP thetaSEXP, SEXP maxorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    rcpp_result_gen = Rcpp::wrap(zpmCP(rho, theta, maxorder));
    return rcpp_result_gen;
END_RCPP
}
// norm_zpm
mat norm_zpm(mat& uzpm, const int& maxorder);
RcppExport SEXP _zernike_norm_zpm(SEXP uzpmSEXP, SEXP maxorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type uzpm(uzpmSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_zpm(uzpm, maxorder));
    return rcpp_result_gen;
END_RCPP
}
// gradzpm
List gradzpm(const vec& x, const vec& y, const int& maxorder, const bool& unit_variance, const bool& return_zpm);
RcppExport SEXP _zernike_gradzpm(SEXP xSEXP, SEXP ySEXP, SEXP maxorderSEXP, SEXP unit_varianceSEXP, SEXP return_zpmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const bool& >::type unit_variance(unit_varianceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type return_zpm(return_zpmSEXP);
    rcpp_result_gen = Rcpp::wrap(gradzpm(x, y, maxorder, unit_variance, return_zpm));
    return rcpp_result_gen;
END_RCPP
}
// zpm_cart
mat zpm_cart(const vec& x, const vec& y, const int& maxorder, const bool& unit_variance);
RcppExport SEXP _zernike_zpm_cart(SEXP xSEXP, SEXP ySEXP, SEXP maxorderSEXP, SEXP unit_varianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type maxorder(maxorderSEXP);
    Rcpp::traits::input_parameter< const bool& >::type unit_variance(unit_varianceSEXP);
    rcpp_result_gen = Rcpp::wrap(zpm_cart(x, y, maxorder, unit_variance));
    return rcpp_result_gen;
END_RCPP
}
// lsfit_qr
vec lsfit_qr(const vec& wf, const mat& zpm);
RcppExport SEXP _zernike_lsfit_qr(SEXP wfSEXP, SEXP zpmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type wf(wfSEXP);
    Rcpp::traits::input_parameter< const mat& >::type zpm(zpmSEXP);
    rcpp_result_gen = Rcpp::wrap(lsfit_qr(wf, zpm));
    return rcpp_result_gen;
END_RCPP
}
// lsfit_norm
vec lsfit_norm(const vec& wf, const mat& zpm);
RcppExport SEXP _zernike_lsfit_norm(SEXP wfSEXP, SEXP zpmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vec& >::type wf(wfSEXP);
    Rcpp::traits::input_parameter< const mat& >::type zpm(zpmSEXP);
    rcpp_result_gen = Rcpp::wrap(lsfit_norm(wf, zpm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_zernike_fitzernikes", (DL_FUNC) &_zernike_fitzernikes, 9},
    {"_zernike_gpcapsiC", (DL_FUNC) &_zernike_gpcapsiC, 4},
    {"_zernike_id_dxy_uw", (DL_FUNC) &_zernike_id_dxy_uw, 7},
    {"_zernike_id_uw", (DL_FUNC) &_zernike_id_uw, 3},
    {"_zernike_lspsiC", (DL_FUNC) &_zernike_lspsiC, 3},
    {"_zernike_aiapsiC", (DL_FUNC) &_zernike_aiapsiC, 5},
    {"_zernike_q_uw", (DL_FUNC) &_zernike_q_uw, 4},
    {"_zernike_readraw", (DL_FUNC) &_zernike_readraw, 2},
    {"_zernike_rescale", (DL_FUNC) &_zernike_rescale, 2},
    {"_zernike_rzernike", (DL_FUNC) &_zernike_rzernike, 3},
    {"_zernike_pwrap", (DL_FUNC) &_zernike_pwrap, 1},
    {"_zernike_tiltpsiC", (DL_FUNC) &_zernike_tiltpsiC, 6},
    {"_zernike_gol_welsch", (DL_FUNC) &_zernike_gol_welsch, 2},
    {"_zernike_rzernike_ann", (DL_FUNC) &_zernike_rzernike_ann, 6},
    {"_zernike_zapm", (DL_FUNC) &_zernike_zapm, 5},
    {"_zernike_zapm_iso", (DL_FUNC) &_zernike_zapm_iso, 5},
    {"_zernike_zapm_128", (DL_FUNC) &_zernike_zapm_128, 5},
    {"_zernike_zapm_iso_128", (DL_FUNC) &_zernike_zapm_iso_128, 5},
    {"_zernike_zpm", (DL_FUNC) &_zernike_zpm, 3},
    {"_zernike_zpmCP", (DL_FUNC) &_zernike_zpmCP, 3},
    {"_zernike_norm_zpm", (DL_FUNC) &_zernike_norm_zpm, 2},
    {"_zernike_gradzpm", (DL_FUNC) &_zernike_gradzpm, 5},
    {"_zernike_zpm_cart", (DL_FUNC) &_zernike_zpm_cart, 4},
    {"_zernike_lsfit_qr", (DL_FUNC) &_zernike_lsfit_qr, 2},
    {"_zernike_lsfit_norm", (DL_FUNC) &_zernike_lsfit_norm, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_zernike(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
