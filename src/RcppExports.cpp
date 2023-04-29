// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eauc
double eauc(arma::vec theta, arma::mat X, arma::vec Y);
RcppExport SEXP _maxAUC_eauc(SEXP thetaSEXP, SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(eauc(theta, X, Y));
    return rcpp_result_gen;
END_RCPP
}
// dqncpp
arma::vec dqncpp(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var);
RcppExport SEXP _maxAUC_dqncpp(SEXP thetaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(dqncpp(theta, X, Y, var));
    return rcpp_result_gen;
END_RCPP
}
// an_anchor
arma::mat an_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var);
RcppExport SEXP _maxAUC_an_anchor(SEXP thetaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(an_anchor(theta, X, Y, var));
    return rcpp_result_gen;
END_RCPP
}
// vn_anchor
arma::mat vn_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var);
RcppExport SEXP _maxAUC_vn_anchor(SEXP thetaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(vn_anchor(theta, X, Y, var));
    return rcpp_result_gen;
END_RCPP
}
// dn_anchor
arma::mat dn_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var);
RcppExport SEXP _maxAUC_dn_anchor(SEXP thetaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(dn_anchor(theta, X, Y, var));
    return rcpp_result_gen;
END_RCPP
}
// newton_raphson_anchor
arma::mat newton_raphson_anchor(arma::vec theta_initial, double tol, int iteration, double gamma, arma::mat X, arma::vec Y, arma::mat var);
RcppExport SEXP _maxAUC_newton_raphson_anchor(SEXP theta_initialSEXP, SEXP tolSEXP, SEXP iterationSEXP, SEXP gammaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta_initial(theta_initialSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(newton_raphson_anchor(theta_initial, tol, iteration, gamma, X, Y, var));
    return rcpp_result_gen;
END_RCPP
}
// eauc_l1
double eauc_l1(arma::vec beta, arma::mat X, arma::vec Y, bool silence);
RcppExport SEXP _maxAUC_eauc_l1(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP silenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type silence(silenceSEXP);
    rcpp_result_gen = Rcpp::wrap(eauc_l1(beta, X, Y, silence));
    return rcpp_result_gen;
END_RCPP
}
// eauc_sort
double eauc_sort(arma::vec beta, arma::mat X, arma::vec Y, bool silence);
RcppExport SEXP _maxAUC_eauc_sort(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP silenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type silence(silenceSEXP);
    rcpp_result_gen = Rcpp::wrap(eauc_sort(beta, X, Y, silence));
    return rcpp_result_gen;
END_RCPP
}
// triang
double triang(double x, double sigma0);
RcppExport SEXP _maxAUC_triang(SEXP xSEXP, SEXP sigma0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    rcpp_result_gen = Rcpp::wrap(triang(x, sigma0));
    return rcpp_result_gen;
END_RCPP
}
// tauc_sort
double tauc_sort(arma::vec beta, arma::mat X, arma::vec Y, double sigma0, bool silence);
RcppExport SEXP _maxAUC_tauc_sort(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP sigma0SEXP, SEXP silenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< bool >::type silence(silenceSEXP);
    rcpp_result_gen = Rcpp::wrap(tauc_sort(beta, X, Y, sigma0, silence));
    return rcpp_result_gen;
END_RCPP
}
// varauc_l1
double varauc_l1(arma::vec beta, arma::mat X, arma::vec Y);
RcppExport SEXP _maxAUC_varauc_l1(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(varauc_l1(beta, X, Y));
    return rcpp_result_gen;
END_RCPP
}
// dtauc_opt
arma::mat dtauc_opt(arma::vec beta, arma::vec beta_k, arma::mat X, arma::vec Y, double sigma0, double w, double t);
RcppExport SEXP _maxAUC_dtauc_opt(SEXP betaSEXP, SEXP beta_kSEXP, SEXP XSEXP, SEXP YSEXP, SEXP sigma0SEXP, SEXP wSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_k(beta_kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type sigma0(sigma0SEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(dtauc_opt(beta, beta_k, X, Y, sigma0, w, t));
    return rcpp_result_gen;
END_RCPP
}
// vn
arma::mat vn(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor);
RcppExport SEXP _maxAUC_vn(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP, SEXP anchorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    Rcpp::traits::input_parameter< int >::type anchor(anchorSEXP);
    rcpp_result_gen = Rcpp::wrap(vn(beta, X, Y, var, anchor));
    return rcpp_result_gen;
END_RCPP
}
// an
arma::mat an(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor);
RcppExport SEXP _maxAUC_an(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP, SEXP anchorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    Rcpp::traits::input_parameter< int >::type anchor(anchorSEXP);
    rcpp_result_gen = Rcpp::wrap(an(beta, X, Y, var, anchor));
    return rcpp_result_gen;
END_RCPP
}
// dn
arma::mat dn(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor);
RcppExport SEXP _maxAUC_dn(SEXP betaSEXP, SEXP XSEXP, SEXP YSEXP, SEXP varSEXP, SEXP anchorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    Rcpp::traits::input_parameter< int >::type anchor(anchorSEXP);
    rcpp_result_gen = Rcpp::wrap(dn(beta, X, Y, var, anchor));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_maxAUC_eauc", (DL_FUNC) &_maxAUC_eauc, 3},
    {"_maxAUC_dqncpp", (DL_FUNC) &_maxAUC_dqncpp, 4},
    {"_maxAUC_an_anchor", (DL_FUNC) &_maxAUC_an_anchor, 4},
    {"_maxAUC_vn_anchor", (DL_FUNC) &_maxAUC_vn_anchor, 4},
    {"_maxAUC_dn_anchor", (DL_FUNC) &_maxAUC_dn_anchor, 4},
    {"_maxAUC_newton_raphson_anchor", (DL_FUNC) &_maxAUC_newton_raphson_anchor, 7},
    {"_maxAUC_eauc_l1", (DL_FUNC) &_maxAUC_eauc_l1, 4},
    {"_maxAUC_eauc_sort", (DL_FUNC) &_maxAUC_eauc_sort, 4},
    {"_maxAUC_triang", (DL_FUNC) &_maxAUC_triang, 2},
    {"_maxAUC_tauc_sort", (DL_FUNC) &_maxAUC_tauc_sort, 5},
    {"_maxAUC_varauc_l1", (DL_FUNC) &_maxAUC_varauc_l1, 3},
    {"_maxAUC_dtauc_opt", (DL_FUNC) &_maxAUC_dtauc_opt, 7},
    {"_maxAUC_vn", (DL_FUNC) &_maxAUC_vn, 5},
    {"_maxAUC_an", (DL_FUNC) &_maxAUC_an, 5},
    {"_maxAUC_dn", (DL_FUNC) &_maxAUC_dn, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_maxAUC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
