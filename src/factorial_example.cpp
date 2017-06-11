#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

//'Calculate the factorial of an integer number
//'
//'Meant to use only for teaching purposes not for replacing the built-in
//'function \code{factorial()}.
//'
//'@param n An integer number
//'
//'@details A C++ implementation using a \code{for} loop.
//'
//'@return The factorial of \code{n}
//'@export
//'
//' @examples
//'    fact_cpp(10);
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

//'Calculate the factorial of an integer number
//'
//'Meant to use only for teaching purposes not for replacing the built-in
//'function \code{factorial()}.
//'
//'@param n An integer number
//'
//'@details A C++ implementation using recursion.
//'
//'@return The factorial of \code{n}
//'@export
//'
//' @examples
//'    fact_cpp(10);
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











