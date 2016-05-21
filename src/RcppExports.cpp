// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getTransitionParameters
DataFrame getTransitionParameters(std::string tpl, std::string mdl, NumericVector snrs);
RcppExport SEXP cc2r_getTransitionParameters(SEXP tplSEXP, SEXP mdlSEXP, SEXP snrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type tpl(tplSEXP);
    Rcpp::traits::input_parameter< std::string >::type mdl(mdlSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type snrs(snrsSEXP);
    __result = Rcpp::wrap(getTransitionParameters(tpl, mdl, snrs));
    return __result;
END_RCPP
}
// getScore
DataFrame getScore(std::string read, std::string tpl, std::string mdl, std::vector<int>& pws, NumericVector snrs);
RcppExport SEXP cc2r_getScore(SEXP readSEXP, SEXP tplSEXP, SEXP mdlSEXP, SEXP pwsSEXP, SEXP snrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type read(readSEXP);
    Rcpp::traits::input_parameter< std::string >::type tpl(tplSEXP);
    Rcpp::traits::input_parameter< std::string >::type mdl(mdlSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type pws(pwsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type snrs(snrsSEXP);
    __result = Rcpp::wrap(getScore(read, tpl, mdl, pws, snrs));
    return __result;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP cc2r_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello());
    return __result;
END_RCPP
}
