// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rms_scan_grid
Rcpp::DataFrame rms_scan_grid(Rcpp::NumericMatrix ref, Rcpp::NumericMatrix mov, Rcpp::NumericMatrix param_grid);
RcppExport SEXP _lidRalignment_rms_scan_grid(SEXP refSEXP, SEXP movSEXP, SEXP param_gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mov(movSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type param_grid(param_gridSEXP);
    rcpp_result_gen = Rcpp::wrap(rms_scan_grid(ref, mov, param_grid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_icp
Rcpp::NumericMatrix cpp_icp(Eigen::MatrixXd& source_mat, Eigen::MatrixXd& target_mat, bool tz_only, bool rz_only, int max_iterations, int overlap, double tolerance);
RcppExport SEXP _lidRalignment_cpp_icp(SEXP source_matSEXP, SEXP target_matSEXP, SEXP tz_onlySEXP, SEXP rz_onlySEXP, SEXP max_iterationsSEXP, SEXP overlapSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type source_mat(source_matSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type target_mat(target_matSEXP);
    Rcpp::traits::input_parameter< bool >::type tz_only(tz_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type rz_only(rz_onlySEXP);
    Rcpp::traits::input_parameter< int >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type overlap(overlapSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_icp(source_mat, target_mat, tz_only, rz_only, max_iterations, overlap, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// cpp_smooth3d
DataFrame cpp_smooth3d(S4 las, NumericVector radius, NumericVector weight, int ncpu, bool pgbar, bool verbose);
RcppExport SEXP _lidRalignment_cpp_smooth3d(SEXP lasSEXP, SEXP radiusSEXP, SEXP weightSEXP, SEXP ncpuSEXP, SEXP pgbarSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type las(lasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type ncpu(ncpuSEXP);
    Rcpp::traits::input_parameter< bool >::type pgbar(pgbarSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_smooth3d(las, radius, weight, ncpu, pgbar, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lidRalignment_rms_scan_grid", (DL_FUNC) &_lidRalignment_rms_scan_grid, 3},
    {"_lidRalignment_cpp_icp", (DL_FUNC) &_lidRalignment_cpp_icp, 7},
    {"_lidRalignment_cpp_smooth3d", (DL_FUNC) &_lidRalignment_cpp_smooth3d, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_lidRalignment(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
