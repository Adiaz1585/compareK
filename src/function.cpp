#include <RcppArmadillo.h>

using namespace Rcpp;

//compareKgroups

List compareKgroups(int ngroup, int nsamples, int burn, IntegerMatrix y, IntegerMatrix missing, int demleader, int repleader, IntegerVector group, int thin);
RcppExport SEXP _pract3_compareKgroups(SEXP ngroupSEXP, SEXP nsamplesSEXP, SEXP burnSEXP, SEXP ySEXP, SEXP missingSEXP, SEXP demleaderSEXP, SEXP repleaderSEXP, SEXP groupSEXP, SEXP thinSEXP)
{
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope Rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ngroup(ngroupSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type missing(missingSEXP);
    Rcpp::traits::input_parameter< int >::type demleader(demleaderSEXP);
    Rcpp::traits::input_parameter< int >::type repleader(repleaderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(compareKgroups(ngroup, nsamples, burn, y, missing, demleader, repleader, group, thin));
    return rcpp_result_gen;
    END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pract3_compareKgroups", (DL_FUNC) _pract3_compareKgroups, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_pract3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
