#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double fact_cpp(double n){
  double f = 1;
  if (n < 0) {
    Rcpp::stop("No factorial for negative numbers");
  } else if (n == 0) {
    return 1;
  } else {
    for (int i = 1; i <= n; ++i) {
      f *= i;
    }
    return f;
  }
}

// [[Rcpp::export]]
double fact_rec_cpp(double n){
  if (n < 0) {
    Rcpp::stop("No factorial for negative numbers");
  } else if (n == 0) {
    return 1;
  } else {
    return n*fact_rec_cpp(n - 1);
  }
}











