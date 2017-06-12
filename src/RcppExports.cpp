// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// gen_fib_cpp
std::vector<double> gen_fib_cpp(const int n);
RcppExport SEXP spsintro_gen_fib_cpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_fib_cpp(n));
    return rcpp_result_gen;
END_RCPP
}
// find_max
double find_max(const std::vector<double>& x);
RcppExport SEXP spsintro_find_max(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_max(x));
    return rcpp_result_gen;
END_RCPP
}
// fact_cpp
double fact_cpp(double n);
RcppExport SEXP spsintro_fact_cpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(fact_cpp(n));
    return rcpp_result_gen;
END_RCPP
}
// fact_rec_cpp
double fact_rec_cpp(double n);
RcppExport SEXP spsintro_fact_rec_cpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(fact_rec_cpp(n));
    return rcpp_result_gen;
END_RCPP
}
// calc_dist
double calc_dist(const Eigen::VectorXd& x_i, const Eigen::VectorXd& x_j, const Eigen::MatrixXd& M);
RcppExport SEXP spsintro_calc_dist(SEXP x_iSEXP, SEXP x_jSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x_i(x_iSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x_j(x_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_dist(x_i, x_j, M));
    return rcpp_result_gen;
END_RCPP
}
// find_knn
Rcpp::List find_knn(const Eigen::MatrixXd& Y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& M, const unsigned& k);
RcppExport SEXP spsintro_find_knn(SEXP YSEXP, SEXP XSEXP, SEXP MSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(find_knn(Y, X, M, k));
    return rcpp_result_gen;
END_RCPP
}
// show_type_ranges
void show_type_ranges();
RcppExport SEXP spsintro_show_type_ranges() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_type_ranges();
    return R_NilValue;
END_RCPP
}
// type_conversions
void type_conversions();
RcppExport SEXP spsintro_type_conversions() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    type_conversions();
    return R_NilValue;
END_RCPP
}
// mix_signed_unsigned
void mix_signed_unsigned();
RcppExport SEXP spsintro_mix_signed_unsigned() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    mix_signed_unsigned();
    return R_NilValue;
END_RCPP
}
// illustrate_scope
void illustrate_scope();
RcppExport SEXP spsintro_illustrate_scope() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    illustrate_scope();
    return R_NilValue;
END_RCPP
}
// show_references
void show_references();
RcppExport SEXP spsintro_show_references() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_references();
    return R_NilValue;
END_RCPP
}
// show_pointers
void show_pointers();
RcppExport SEXP spsintro_show_pointers() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    show_pointers();
    return R_NilValue;
END_RCPP
}
// create_tree_struct
void create_tree_struct(double dbh, double height, int no_stems, std::string spec_name);
RcppExport SEXP spsintro_create_tree_struct(SEXP dbhSEXP, SEXP heightSEXP, SEXP no_stemsSEXP, SEXP spec_nameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< double >::type height(heightSEXP);
    Rcpp::traits::input_parameter< int >::type no_stems(no_stemsSEXP);
    Rcpp::traits::input_parameter< std::string >::type spec_name(spec_nameSEXP);
    create_tree_struct(dbh, height, no_stems, spec_name);
    return R_NilValue;
END_RCPP
}
// to_upper_range_for
std::string to_upper_range_for(std::string str);
RcppExport SEXP spsintro_to_upper_range_for(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type str(strSEXP);
    rcpp_result_gen = Rcpp::wrap(to_upper_range_for(str));
    return rcpp_result_gen;
END_RCPP
}
// to_upper_subscript
std::string to_upper_subscript(std::string str);
RcppExport SEXP spsintro_to_upper_subscript(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type str(strSEXP);
    rcpp_result_gen = Rcpp::wrap(to_upper_subscript(str));
    return rcpp_result_gen;
END_RCPP
}
// add_elem
std::vector<int> add_elem(int n);
RcppExport SEXP spsintro_add_elem(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(add_elem(n));
    return rcpp_result_gen;
END_RCPP
}
// sum_elem
double sum_elem(std::vector<double> v);
RcppExport SEXP spsintro_sum_elem(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_elem(v));
    return rcpp_result_gen;
END_RCPP
}
// squ_elem
std::vector<double> squ_elem(std::vector<double> v);
RcppExport SEXP spsintro_squ_elem(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(squ_elem(v));
    return rcpp_result_gen;
END_RCPP
}
// squ_elem_it
std::vector<double> squ_elem_it(std::vector<double> v);
RcppExport SEXP spsintro_squ_elem_it(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(squ_elem_it(v));
    return rcpp_result_gen;
END_RCPP
}
// arr_access
void arr_access();
RcppExport SEXP spsintro_arr_access() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    arr_access();
    return R_NilValue;
END_RCPP
}
// ret_rem
int ret_rem(int a, int b);
RcppExport SEXP spsintro_ret_rem(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(ret_rem(a, b));
    return rcpp_result_gen;
END_RCPP
}
// dem_cond
void dem_cond(int x);
RcppExport SEXP spsintro_dem_cond(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    dem_cond(x);
    return R_NilValue;
END_RCPP
}
// dem_if_else
int dem_if_else(double dbh);
RcppExport SEXP spsintro_dem_if_else(SEXP dbhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type dbh(dbhSEXP);
    rcpp_result_gen = Rcpp::wrap(dem_if_else(dbh));
    return rcpp_result_gen;
END_RCPP
}
// find_first_neg
int find_first_neg(std::vector<double> x);
RcppExport SEXP spsintro_find_first_neg(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_first_neg(x));
    return rcpp_result_gen;
END_RCPP
}
// find_first_neg_for
int find_first_neg_for(std::vector<double> x);
RcppExport SEXP spsintro_find_first_neg_for(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_first_neg_for(x));
    return rcpp_result_gen;
END_RCPP
}
// find_first_neg_rf
int find_first_neg_rf(std::vector<double> x);
RcppExport SEXP spsintro_find_first_neg_rf(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_first_neg_rf(x));
    return rcpp_result_gen;
END_RCPP
}
// find_first_neg_dw
unsigned find_first_neg_dw(std::vector<double> x);
RcppExport SEXP spsintro_find_first_neg_dw(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(find_first_neg_dw(x));
    return rcpp_result_gen;
END_RCPP
}
// round_cpp
double round_cpp(double x);
RcppExport SEXP spsintro_round_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(round_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// use_polygon_class
void use_polygon_class();
RcppExport SEXP spsintro_use_polygon_class() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    use_polygon_class();
    return R_NilValue;
END_RCPP
}
// inheritance
void inheritance();
RcppExport SEXP spsintro_inheritance() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    inheritance();
    return R_NilValue;
END_RCPP
}
// bubble_sort
std::vector<double> bubble_sort(std::vector<double>& x);
RcppExport SEXP spsintro_bubble_sort(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double>& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bubble_sort(x));
    return rcpp_result_gen;
END_RCPP
}
// use_matrix
Rcpp::NumericMatrix use_matrix(const Rcpp::NumericMatrix& x);
RcppExport SEXP spsintro_use_matrix(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(use_matrix(x));
    return rcpp_result_gen;
END_RCPP
}
// use_list
Rcpp::List use_list(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& y);
RcppExport SEXP spsintro_use_list(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(use_list(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calc_poly_area
double calc_poly_area(const std::vector<double>& x, const std::vector<double>& y);
RcppExport SEXP spsintro_calc_poly_area(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(calc_poly_area(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"spsintro_gen_fib_cpp", (DL_FUNC) &spsintro_gen_fib_cpp, 1},
    {"spsintro_find_max", (DL_FUNC) &spsintro_find_max, 1},
    {"spsintro_fact_cpp", (DL_FUNC) &spsintro_fact_cpp, 1},
    {"spsintro_fact_rec_cpp", (DL_FUNC) &spsintro_fact_rec_cpp, 1},
    {"spsintro_calc_dist", (DL_FUNC) &spsintro_calc_dist, 3},
    {"spsintro_find_knn", (DL_FUNC) &spsintro_find_knn, 4},
    {"spsintro_show_type_ranges", (DL_FUNC) &spsintro_show_type_ranges, 0},
    {"spsintro_type_conversions", (DL_FUNC) &spsintro_type_conversions, 0},
    {"spsintro_mix_signed_unsigned", (DL_FUNC) &spsintro_mix_signed_unsigned, 0},
    {"spsintro_illustrate_scope", (DL_FUNC) &spsintro_illustrate_scope, 0},
    {"spsintro_show_references", (DL_FUNC) &spsintro_show_references, 0},
    {"spsintro_show_pointers", (DL_FUNC) &spsintro_show_pointers, 0},
    {"spsintro_create_tree_struct", (DL_FUNC) &spsintro_create_tree_struct, 4},
    {"spsintro_to_upper_range_for", (DL_FUNC) &spsintro_to_upper_range_for, 1},
    {"spsintro_to_upper_subscript", (DL_FUNC) &spsintro_to_upper_subscript, 1},
    {"spsintro_add_elem", (DL_FUNC) &spsintro_add_elem, 1},
    {"spsintro_sum_elem", (DL_FUNC) &spsintro_sum_elem, 1},
    {"spsintro_squ_elem", (DL_FUNC) &spsintro_squ_elem, 1},
    {"spsintro_squ_elem_it", (DL_FUNC) &spsintro_squ_elem_it, 1},
    {"spsintro_arr_access", (DL_FUNC) &spsintro_arr_access, 0},
    {"spsintro_ret_rem", (DL_FUNC) &spsintro_ret_rem, 2},
    {"spsintro_dem_cond", (DL_FUNC) &spsintro_dem_cond, 1},
    {"spsintro_dem_if_else", (DL_FUNC) &spsintro_dem_if_else, 1},
    {"spsintro_find_first_neg", (DL_FUNC) &spsintro_find_first_neg, 1},
    {"spsintro_find_first_neg_for", (DL_FUNC) &spsintro_find_first_neg_for, 1},
    {"spsintro_find_first_neg_rf", (DL_FUNC) &spsintro_find_first_neg_rf, 1},
    {"spsintro_find_first_neg_dw", (DL_FUNC) &spsintro_find_first_neg_dw, 1},
    {"spsintro_round_cpp", (DL_FUNC) &spsintro_round_cpp, 1},
    {"spsintro_use_polygon_class", (DL_FUNC) &spsintro_use_polygon_class, 0},
    {"spsintro_inheritance", (DL_FUNC) &spsintro_inheritance, 0},
    {"spsintro_bubble_sort", (DL_FUNC) &spsintro_bubble_sort, 1},
    {"spsintro_use_matrix", (DL_FUNC) &spsintro_use_matrix, 1},
    {"spsintro_use_list", (DL_FUNC) &spsintro_use_list, 2},
    {"spsintro_calc_poly_area", (DL_FUNC) &spsintro_calc_poly_area, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spsintro(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
