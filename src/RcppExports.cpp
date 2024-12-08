// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// distanceMatrixSmoothing
NumericMatrix distanceMatrixSmoothing(NumericMatrix dmat, NumericMatrix umat);
RcppExport SEXP _ecosmooth_distanceMatrixSmoothing(SEXP dmatSEXP, SEXP umatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dmat(dmatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type umat(umatSEXP);
    rcpp_result_gen = Rcpp::wrap(distanceMatrixSmoothing(dmat, umat));
    return rcpp_result_gen;
END_RCPP
}
// rectangularMatrixSmoothing
NumericMatrix rectangularMatrixSmoothing(NumericMatrix smat, NumericMatrix umat);
RcppExport SEXP _ecosmooth_rectangularMatrixSmoothing(SEXP smatSEXP, SEXP umatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type smat(smatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type umat(umatSEXP);
    rcpp_result_gen = Rcpp::wrap(rectangularMatrixSmoothing(smat, umat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ecosmooth_distanceMatrixSmoothing", (DL_FUNC) &_ecosmooth_distanceMatrixSmoothing, 2},
    {"_ecosmooth_rectangularMatrixSmoothing", (DL_FUNC) &_ecosmooth_rectangularMatrixSmoothing, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ecosmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
