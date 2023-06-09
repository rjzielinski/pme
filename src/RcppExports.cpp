// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// norm_euclidean
double norm_euclidean(arma::vec x);
RcppExport SEXP _pme_norm_euclidean(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_euclidean(x));
    return rcpp_result_gen;
END_RCPP
}
// dist_euclidean
double dist_euclidean(arma::vec x, arma::vec y);
RcppExport SEXP _pme_dist_euclidean(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_euclidean(x, y));
    return rcpp_result_gen;
END_RCPP
}
// eta_kernel
double eta_kernel(arma::vec t, int lambda);
RcppExport SEXP _pme_eta_kernel(SEXP tSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_kernel(t, lambda));
    return rcpp_result_gen;
END_RCPP
}
// calcE
arma::mat calcE(arma::mat x, int lambda);
RcppExport SEXP _pme_calcE(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(calcE(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// etaFunc
arma::mat etaFunc(arma::vec t, arma::mat tau, int lambda);
RcppExport SEXP _pme_etaFunc(SEXP tSEXP, SEXP tauSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(etaFunc(t, tau, lambda));
    return rcpp_result_gen;
END_RCPP
}
// calc_nearest_x
arma::vec calc_nearest_x(arma::mat df, arma::mat x);
RcppExport SEXP _pme_calc_nearest_x(SEXP dfSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_nearest_x(df, x));
    return rcpp_result_gen;
END_RCPP
}
// calc_init_param
arma::mat calc_init_param(arma::mat df, arma::mat tnew, arma::vec nearest_x);
RcppExport SEXP _pme_calc_init_param(SEXP dfSEXP, SEXP tnewSEXP, SEXP nearest_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tnew(tnewSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nearest_x(nearest_xSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_init_param(df, tnew, nearest_x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pme_norm_euclidean", (DL_FUNC) &_pme_norm_euclidean, 1},
    {"_pme_dist_euclidean", (DL_FUNC) &_pme_dist_euclidean, 2},
    {"_pme_eta_kernel", (DL_FUNC) &_pme_eta_kernel, 2},
    {"_pme_calcE", (DL_FUNC) &_pme_calcE, 2},
    {"_pme_etaFunc", (DL_FUNC) &_pme_etaFunc, 3},
    {"_pme_calc_nearest_x", (DL_FUNC) &_pme_calc_nearest_x, 2},
    {"_pme_calc_init_param", (DL_FUNC) &_pme_calc_init_param, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_pme(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
