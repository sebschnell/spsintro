#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix use_matrix(const Rcpp::NumericMatrix &x) {
  Rcpp::NumericMatrix m1(x);
  Rcpp::NumericMatrix m2(5, 4);

  Rcpp::Rcout << m1(1, 2) << std::endl;

  return m2;
}

// [[Rcpp::export]]
Rcpp::List use_list(const Rcpp::NumericMatrix &x, const Rcpp::NumericVector &y) {
  return Rcpp::List::create(Rcpp::Named("m") = x,
                            Rcpp::Named("v") = y);
}
